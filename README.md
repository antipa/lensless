# lensless
Reconstruct images from lensless coding cameras

## 2D deconvolution using proximal gradient descent
The script `realtime_deconv.m` takes care of loading in measured data, preparing file names etc., then calls `diffuser_2d_deconv_v2.m` to do the deconvolution. `diffuser_2d_deconv_v2` relies on the code in `antipa/proxMin` to function, so you'll need that. At the top of `realtime_deconv` are:
1. `input_folder`  where you'll need to specify the path to the image to process 
2. `camera_type` string containing camera type. Currently needs to be either `'pco'` or `'flea3'`. This must be set because each camera's raw data contains a different bias, which, if not properly subtracted, leads to poor reconstruction.
3. `process_color` String specifying which color channel to use. For `flea3`, mono is the only option, but for `pco`, choices are:
  * `'mono'`: average RGB data in PSF and raw measurement
  * `'red'`: use only red channel (after demosaicing)
  * `'green'`: green only
  * `'blue'`: blue
  * `'all'`: solve all 3 colors independently, then create RGB output

At the top of `diffuser_2d_deconv_v2`, you can specify `ds` to be an integer. The images will all be downsampled (using a 'box' antialias filter) by this amount.

You can also set the regularizer parameters in the code, but this needs to be made more user friendly before it is recommended.
