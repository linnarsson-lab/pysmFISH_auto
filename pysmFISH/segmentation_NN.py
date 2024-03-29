import numpy as np
import matplotlib.pyplot as plt


class Segmenation_NN:
    def __init__(self, mode: str = "stardist", diameter: float = 25, model=None):
        """Instantiate NeuralNet segmentation class.

        Args:
            mode (str, optional): Segmentation algorithm to use. Choose from
                "stardist" or "cellpose", Defaults to 'stardist'.
            diameter (int, optional): Only needed for Cellpose. Cell diameter
                prior in pixels. If None is given, Cellpose will calculate the
                diameter. Defaults to 25, which is a value that worked for mouse
                adult brain nuclei imaged with the 40X objective.
                Defaults to 25.
        """
        self.mode = mode
        self.model = model
        self.diameter = diameter

        if self.mode == "stardist":
            from stardist.models import StarDist2D, Config2D

            self.StarDist2D = StarDist2D
            self.Config2D = Config2D
            from csbdeep.utils import normalize

            self.stardist_normalize = normalize

        if self.mode == "cellpose":
            from cellpose import models

            self.cellpose_models = models
            print(models)

        self.diameter = diameter

    # def cellpose_init_model(
    #     self,
    #     gpu: bool = False,
    #     model_type: str = "nuclei",
    #     net_avg: bool = True,
    #     device: object = None,
    # ) -> None:

    #     """Initiate the Cellpose model.

    #     Initiates the model and returns it.
    #     Refer to Cellpose documentation for more details on input.

    #     https://doi.org/10.1038/s41592-020-01018-x

    #     Args:
    #         gpu (bool, optional): True if Cuda is installed and is to be used.
    #             Defaults to False.
    #         model_type (str, optional): "nuclei" or "cytoplasm" segmentation.
    #             Defaults to 'nuclei'.
    #         net_avg (bool, optional): Averages build-in networks. See Cellpose
    #             documentation. Defaults to True.

    #         device (object, optional): Use saved model. See Cellpose
    #             documentation. Defaults to None.
    #     """
    #     model = self.cellpose_models.Cellpose(
    #         gpu=gpu, model_type=model_type, net_avg=net_avg, device=device
    #     )
    #     return model

    def segment_cellpose(self, image, gpu=False, model_type="nuclei"):
        """Segment image with CellPose.

        Saves segmentation masks as numpy arrays in the output folder. Other
        output of cellpose.model.eval is not saved.

        Args:
            image (np.ndarray): Numpy array with the image of the nuclei.
            image_id (str): String identifier of the image. Image will be saved
                as <image_id>_segmentation_mask.npy.
            output_folder (str): Path to output folder.

            diameter (float, optional): Cell diameter prior in pixels. If None
                is given, Cellpose will calculate the diameter. Defaults to 25,
                which is a value that worked for mouse adult brain nuclei imaged
                with the 40X objective. Defaults to 25.

            model (cellpose model): Instantiated CellPose model. If using this

                function in parallel make a deepcopy of the model.
            gpu (bool, optional): True if Cuda is installed and is to be used
                and if model is not provided directly. Defaults to False.
            model_type (str, optional): "nuclei" or "cytoplasm" segmentation.
                Only if model is not provided directly.
                Defaults to 'nuclei'.


        """
        if self.model == None:
            self.model = self.cellpose_models.Cellpose(
                gpu=gpu, model_type=model_type, net_avg=None, device=True
            )
        # Segment
        mask, flow, style, diam = self.model.eval(
            image, diameter=self.diameter, channels=[0, 0]
        )

        # Return
        return mask

    def segment_stardist(self, image):
        """Segment image with Stardist.

        Saves segmentation masks as numpy arrays in the output folder.


        Args:
            image (np.ndarray): Numpy array with the image of the nuclei.
            model (str): Stardist model to use.
            diameter: Placeholder.


        """
        # Instantiate model
        if self.model == None:
            #Explicitly use CPU
            conf = self.Config2D(use_gpu=False)
            self.model = self.StarDist2D(config=conf).from_pretrained("2D_versatile_fluo")

        # Segment test
        img = self.stardist_normalize(image)
        if isinstance(img, np.ndarray):
            mask, _ = self.model.predict_instances(img)
            return mask
        else:
            return img
        # mask, _ = self.model.predict_instances(self.stardist_normalize(image))

        # Retrun
        return mask

    def segment(self, image):
        """Run segmentation with chosen method

        Args:
            image (np.ndarray): Numpy array with the image of the nuclei.

        Returns:
            np.ndarray : Resulting label mask
        """

        if self.mode == "cellpose":
            return self.segment_cellpose(image, self.diameter)

        if self.mode == "stardist":
            return self.segment_stardist(image)

    def segmentation_inspect(
        self,
        image: np.ndarray,
        mask: np.ndarray,
        figsize: float = 15,
        vmin: float = None,
        vmax: float = None,
        alpha: float = 0.5,
        save=False,
        save_name="segmenation_inspection.png",
        dpi=300,
    ):

        """Plot segmentation mask over input image.

        Args:
            image (np.ndarray): Original image.
            mask (np.ndarray): Segmentation mask.
            figsize (float, optional): Figure size. Figures will be square.
                Defaults to 15.
            vmin (float, optional): Vmin of the image. Image is normalized
                between 0-1. Defaults to None.

            vmax (float, optional): Vmax of the image. Image is normalized

                between 0-1. Defaults to None.
            alpha (float, optional): Alpha value of the mask. Defaults to 0.5.
        """

        # Shuffle the labels to help visualization

        unique = np.unique(mask)
        unique = np.delete(unique, np.array([0]))
        shuffeled = unique.copy()
        np.random.shuffle(shuffeled)
        mix = dict(zip(unique, shuffeled))
        mix[0] = 0

        shuffeled_mask = np.array([mix[i] for i in mask.ravel()])
        shuffeled_mask = shuffeled_mask.reshape(mask.shape)
        shuffeled_mask = np.ma.masked_where(shuffeled_mask == 0, shuffeled_mask)

        # Rescale image
        image = image / image.max()

        # Plot
        fig = plt.figure(figsize=(figsize, figsize))
        plt.imshow(image, cmap=plt.cm.gray, vmin=vmin, vmax=vmax)
        plt.imshow(shuffeled_mask, cmap="hsv", alpha=alpha)
        if save:
            plt.savefig(save_name, dpi=dpi)
