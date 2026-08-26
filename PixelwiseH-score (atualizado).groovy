import groovy.transform.CompileStatic
import org.bytedeco.opencv.global.opencv_core
import org.bytedeco.opencv.opencv_core.Mat
import qupath.lib.images.servers.PixelCalibration
import qupath.lib.objects.classes.PathClass
import qupath.lib.regions.ImageRegion
import qupath.opencv.ml.pixel.PixelClassifierTools
import qupath.opencv.ml.pixel.PixelClassifiers
import qupath.opencv.ops.ImageOp
import qupath.opencv.ops.ImageOps
import qupath.opencv.tools.OpenCVTools

import static qupath.lib.gui.scripting.QPEx.*


// SET THESE
DABthresholds = [0.05, 0.15, 0.5, 1] // 0+, 1+, 2+, 3+
Hthreshold = 0.08 // Nuclei in Hematoxylin
double scale = 1 // bigger = downsample for speed

// Pre-request tiles for speed
boolean prerequestTiles = true

// Handle multiple annotations in parallel
boolean parallelizeAnnotations = false


// Make pixel_classifiers folder
String pixelClassifierFolder = buildFilePath(
        PROJECT_BASE_DIR,
        'classifiers',
        'pixel_classifiers'
)
mkdirs(pixelClassifierFolder)


// Check we have objects
def annotations = new ArrayList<>(getSelectedObjects())

if (annotations.isEmpty()) {
    println "Please select some objects!"
    return
}


// Get pixel calibration and current stain vectors
def imageData = getCurrentImageData()
PixelCalibration cal = imageData.getServer().getPixelCalibration()

def stains = getCurrentImageData().getColorDeconvolutionStains()


// Create an ImageDataOp to preprocess and apply thresholds
def op = ImageOps.buildImageDataOp().appendOps(

        ImageOps.Channels.deconvolve(stains),

        ImageOps.Channels.extract(0, 1), // H=0, DAB=1

        new HScoreThresholdOp()
                .zeroPlusThreshold(DABthresholds[0])
                .lowThreshold(DABthresholds[1])
                .mediumThreshold(DABthresholds[2])
                .highThreshold(DABthresholds[3])
                .minStainThreshold(Hthreshold)
)


// Create classes
def classmap = [
        255: PathClass.getInstance("Ignore*"),
        0: PathClass.getInstance("Negative"),
        1: PathClass.getInstance("0+"),
        2: PathClass.getInstance("1+"),
        3: PathClass.getInstance("2+"),
        4: PathClass.getInstance("3+")
]


// Create the classifier
def pixelClassifier = PixelClassifiers.createClassifier(
        op,
        cal.createScaledInstance(scale, scale),
        classmap
)


// Create a pixel classifier server and manager
def pixelServer =
        PixelClassifierTools.createPixelClassificationServer(
                imageData,
                pixelClassifier
        )

def manager =
        PixelClassifierTools.createMeasurementManager(pixelServer)

String measurementPrefix = "H-score"


// Measure the area of each class for each annotation
def annotationStream = annotations.stream()

if (parallelizeAnnotations)
    annotationStream = annotationStream.parallel()

annotationStream.forEach { annotation ->

    if (prerequestTiles) {

        def roi = annotation.getROI()
        def region = ImageRegion.createInstance(roi)

        def tiles = pixelServer
                .getTileRequestManager()
                .getAllTileRequests()
                .findAll(tile ->
                        tile.getRegionRequest().intersects(region)
                )

        tiles.parallelStream().forEach { tile ->
            pixelServer.readRegion(tile.getRegionRequest())
        }
    }

    PixelClassifierTools.addMeasurements(
            [annotation],
            manager,
            measurementPrefix
    )
}


// Calculate H-score
for (def annotation in annotations) {

    double area0plus =
            annotation.measurements[
                    "$measurementPrefix: 0+ area µm^2"
            ]

    double area1 =
            annotation.measurements[
                    "$measurementPrefix: 1+ area µm^2"
            ]

    double area2 =
            annotation.measurements[
                    "$measurementPrefix: 2+ area µm^2"
            ]

    double area3 =
            annotation.measurements[
                    "$measurementPrefix: 3+ area µm^2"
            ]

    double areaNegative =
            annotation.measurements[
                    "$measurementPrefix: Negative area µm^2"
            ]

    double areaDenom =
            area0plus +
            area1 +
            area2 +
            area3 +
            areaNegative

    annotation.measurements['Pixelwise H-score'] =
            (area1 * 1 +
             area2 * 2 +
             area3 * 3) /
            areaDenom * 100
}


/**
 * Custom op to help with computing a pixelwise H-score
 */
class HScoreThresholdOp implements ImageOp {

    double minThreshold = Double.NEGATIVE_INFINITY

    double t0 = 0.05
    double t1 = 0.15
    double t2 = 0.5
    double t3 = 1


    HScoreThresholdOp minStainThreshold(double minThreshold) {
        this.minThreshold = minThreshold
        return this
    }


    // Threshold separating 0 from 0+
    HScoreThresholdOp zeroPlusThreshold(double t0) {
        this.t0 = t0
        return this
    }


    // Threshold separating 0+ from 1+
    HScoreThresholdOp lowThreshold(double t1) {
        this.t1 = t1
        return this
    }


    // Threshold separating 1+ from 2+
    HScoreThresholdOp mediumThreshold(double t2) {
        this.t2 = t2
        return this
    }


    // Threshold separating 2+ from 3+
    HScoreThresholdOp highThreshold(double t3) {
        this.t3 = t3
        return this
    }


    public Mat apply(Mat input) {

        def split = OpenCVTools.splitChannels(input)

        // Apply thresholds to DAB channel
        def matDAB = split[1]

        OpenCVTools.apply(
                matDAB,
                this::applyDABThresholds
        )


        // Apply threshold to hematoxylin channel
        def matH = split[0]

        OpenCVTools.apply(
                matH,
                this::applyHematoxylinThreshold
        )


        // If 0 in both DAB and hematoxylin channels,
        // set to 255 (exclude)
        def matMask = opencv_core.equals(
                opencv_core.max(matH, matDAB).asMat(),
                0.0
        ).asMat()

        matDAB.setTo(
                OpenCVTools.scalarMat(
                        255.0,
                        matDAB.type()
                ),
                matMask
        )

        return matDAB
    }


    @CompileStatic
    public double applyHematoxylinThreshold(double value) {

        if (value < minThreshold)
            return 0

        return 1
    }


    @CompileStatic
    public double applyDABThresholds(double value) {

        if (value < t0)
            return 0

        else if (value < t1)
            return 1

        else if (value < t2)
            return 2

        else if (value < t3)
            return 3

        else
            return 4
    }
}
