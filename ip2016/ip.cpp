#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector>



/*
 * convolve with a box filter
 */
Image* ip_blur_box (Image* src, int size)
{
    double **kernel = new double*[size];
    for (int x=0; x<size; x++){
        kernel[x] = new double[size];
    }
    for(int i=0;i<size;i++){
        for(int j=0; j<size; j++){
            kernel[i][j] = 1.0/(size*size);
        }
    }
    return ip_convolve(src, size, kernel);
}

/*
 * convolve with a gaussian filter
 */
Image* ip_blur_gaussian (Image* src, int size, double sigma)
{
    double **kernel = new double*[size];
    for (int x=0; x<size; x++){
        kernel[x] = new double[size];
    }
    double sum = 0.0;
    for(int i=0;i<size;i++){
        for(int j=0; j<size; j++){
            int span = size / 2;
            int a = i - span;
            int b = j - span;
            double pi = atan(1.0)*4;
            double n = -(a*a + b*b) / (2.0 * sigma * sigma);
            double g_value = (1 / sqrt(2.0 * pi * sigma * sigma)) * exp(n);
            kernel[i][j] = g_value;
            sum += kernel[i][j];
        }
    }
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            kernel[i][j] /= sum;
        }
    }
    return ip_convolve(src, size, kernel);
}


/*
 * convolve with a triangle filter
 */
Image* ip_blur_triangle (Image* src, int size)
{
	
    cerr << "This filter is no longer in the assignment.\n";
    return NULL;
}


/*
 * interpolate with a black image
 */
Image* ip_brighten (Image* src, double alpha)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* black = new Image(width,height);
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            black->setPixel_(w, h, 0, 0);
            black->setPixel_(w, h, 1, 0);
            black->setPixel_(w, h, 2, 0);
        }
    }
    dest = ip_interpolate (src, black, alpha);
    return dest;
}

Image* ip_color_shift(Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            dest->setPixel(w, h, 0, b);
            dest->setPixel(w, h, 1, r);
            dest->setPixel(w, h, 2, g);
        }
    }
    return dest;
}

/*
 * interpolate with the average intensity of the src image
 */
Image* ip_contrast (Image* src, double alpha)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* grey = new Image(width,height);
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            grey->setPixel_(w, h, 0, 0.5);
            grey->setPixel_(w, h, 1, 0.5);
            grey->setPixel_(w, h, 2, 0.5);
        }
    }
    dest = ip_interpolate (src, grey, alpha);
    return dest;
}


/*
 * convolve an image with another image
 */
Image* ip_convolve (Image* src, int size, double** kernel)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image(width,height);
    int hspan = (size-1)/2;
    int wspan = (size-1)/2;
//    cout << "kernel height: " << sizeof(kernel[0]) << endl;
//    cout << "kernel width: " << sizeof(kernel) << endl;
    for (int w=0; w<width;w++){
        for (int h=0; h<height; h++){
            double r = 0; double g = 0; double b = 0;
            for(int dw=-wspan; dw<wspan; dw++){
                for(int dh=-hspan; dh<hspan; dh++){
                    int ww = w + dw; int hh = h + dh;
                    if(ww>=0 && ww<width & hh>=0 & hh<height){
//                        cout << "ww-hh" << ww << " - " << hh << endl;
                        double weight = kernel[dw+wspan][dh+hspan];
//                        cout << "kernel: " << dw+wspan << " - " << dh+hspan << endl;
//                        cout << "weight: " << weight << endl;
                        double rr = src->getPixel(ww,hh,0);
                        double gg = src->getPixel(ww,hh,1);
                        double bb = src->getPixel(ww,hh,2);
                        r = r + weight * rr;
                        g = g + weight * gg;
                        b = b + weight * bb;
                    }
                }
            }
            output->setPixel_(w,h,0,r);
            output->setPixel_(w,h,1,g);
            output->setPixel_(w,h,2,b);
//            cout << "w-h" << w << " - " << h << endl;
        }
    }
    return output;
}


/*
 * use a mask image for a per-pixel alpha value to perform
 * interpolation with a second image
 */
Image* ip_composite(Image* src1, Image* src2, Image* mask)
{
    int width = src1->getWidth();
    int height = src1->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r1 = src1->getPixel(w, h, 0);
            double g1 = src1->getPixel(w, h, 1);
            double b1 = src1->getPixel(w, h, 2);
            double r2 = src2->getPixel(w, h, 0);
            double g2 = src2->getPixel(w, h, 1);
            double b2 = src2->getPixel(w, h, 2);
            double alpha_r = mask->getPixel(w, h, 0);
            double alpha_g = mask->getPixel(w, h, 1);
            double alpha_b = mask->getPixel(w, h, 2);
            double rd = alpha_r*r1 + (1-alpha_r)*r2;
            double gd = alpha_g*g1 + (1-alpha_g)*g2;
            double bd = alpha_b*b1 + (1-alpha_b)*b2;
            dest->setPixel_(w, h, 0, rd);
            dest->setPixel_(w, h, 1, gd);
            dest->setPixel_(w, h, 2, bd);
        }
    }
    return dest;
}


/*
 * cut away all but a subset of the image
 */
Image* ip_crop (Image* src, int x0, int y0, int x1, int y1)
{
    cerr << "This filter has not been implemented 2.\n";
    return NULL;
}


/*
 * convolve with an edge detection kernel
 */
Image* ip_edge_detect (Image* src)
{
    double **kernel = new double*[3];
    for (int x=0; x<3; x++){
        kernel[x] = new double[3];
    }
    for(int i=0;i<3;i++){
        for(int j=0; j<3; j++){
            kernel[i][j] = (i == 1 && j == 1) ? 8 : -1;
        }
    }
    return ip_convolve(src, 3, kernel);
}


/*
 * create a new image with color values from one channel of src
 */
Image* ip_extract (Image* src, int channel)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0; h<height; h++){
            for (int c=0; c<3; c++) {
                if (c != channel) {
                    dest->setPixel(w, h, c, 0);
                } else {
                    double cvalue = src->getPixel(w, h, channel);
                    dest->setPixel(w, h, c, cvalue);
                }
            }
        }
    }
    return dest;
}

/*
 * perform do what you like
 * you must query user for parameters here
 */
Image* ip_fun_warp (Image* src,int samplingMethod)
{
    cerr << "This filter has not been implemented 3.\n";
    return NULL;
}


/*
 * create a new image with values equal to the psychosomatic intensities
 * of the source image
 */
Image* ip_grey (Image* src)
{
    double cr=  0.2126;
    double cg=  0.7152;
    double cb=  0.0722;
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            double gvalue = cr*r + cg*g + cb*b;
            dest->setPixel(w, h, 0, gvalue);
            dest->setPixel(w, h, 1, gvalue);
            dest->setPixel(w, h, 2, gvalue);
        }
    }
    return dest;
}

/*
 * shift image by dx, dy
 *
*/
Image* ip_image_shift(Image* src, int dx, int dy)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            int nw = (w+dx) % width;
            int nh = (h+dy) % height;
            double r = src->getPixel(nw, nh, 0);
            double g = src->getPixel(nw, nh, 1);
            double b = src->getPixel(nw, nh, 2);
            dest->setPixel(w, h, 0, r);
            dest->setPixel(w, h, 1, g);
            dest->setPixel(w, h, 2, b);
        }
    }
    return dest;
}

/*
 * interpolate an image with another image
 */
Image* ip_interpolate (Image* src1, Image* src2, double alpha)
{
    int width = src1->getWidth();
    int height = src1->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r1 = src1->getPixel(w, h, 0);
            double g1 = src1->getPixel(w, h, 1);
            double b1 = src1->getPixel(w, h, 2);
            double r2 = src2->getPixel(w, h, 0);
            double g2 = src2->getPixel(w, h, 1);
            double b2 = src2->getPixel(w, h, 2);
            double rd = alpha*r1 + (1-alpha)*r2;
            double gd = alpha*g1 + (1-alpha)*g2;
            double bd = alpha*b1 + (1-alpha)*b2;
            dest->setPixel_(w, h, 0, rd);
            dest->setPixel_(w, h, 1, gd);
            dest->setPixel_(w, h, 2, bd);
        }
    }
    return dest;
}

/*
 * subtract the image from a white image
 */
Image* ip_invert (Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    Image* grey = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            grey->setPixel_(w, h, 0, 0.5);
            grey->setPixel_(w, h, 1, 0.5);
            grey->setPixel_(w, h, 2, 0.5);
        }
    }
    dest = ip_interpolate (src, grey, -1);
    return dest;
}

/*
 * median filter
 */
Image* ip_median (Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width-2; w++){
        for (int h=0;h<height-2; h++){
            // construct vector for sorting
            vector<double> vec_r;
            vector<double> vec_g;
            vector<double> vec_b;
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    // get r,g,b values
                    double r = src->getPixel(w+i, h+j, 0);
                    double g = src->getPixel(w+i, h+j, 1);
                    double b = src->getPixel(w+i, h+j, 2);
                    vec_r.push_back(r);
                    vec_g.push_back(g);
                    vec_b.push_back(b);
                }
            }
            // sorting r,g,b vectors
            sort(vec_r.begin(), vec_r.end());
            sort(vec_g.begin(), vec_g.end());
            sort(vec_b.begin(), vec_b.end());
            dest->setPixel(w, h, 0, vec_r[4]);
            dest->setPixel(w, h, 1, vec_g[4]);
            dest->setPixel(w, h, 2, vec_b[4]);
        }
    }
    return dest;
}

/* misc
*/
Image* ip_misc(Image* src) 
{
    cerr << "This filter has not been implemented 4.\n";
    return NULL;
}


/*
 * round each pixel to the nearest value in the new number of bits
 */
Image* ip_quantize_simple (Image* src, int bitsPerChannel)
{
    cerr << "This filter has not been implemented 5.\n";
    return NULL;
}


/*
 * dither each pixel to the nearest value in the new number of bits
 * using a static 2x2 matrix
 */
Image* ip_quantize_ordered (Image* src, int bitsPerChannel)
{
    cerr << "This filter has not been implemented 6.\n";
    return NULL;
}


/*
 * dither each pixel to the nearest value in the new number of bits
 * using error diffusion
 */
Image* ip_quantize_fs (Image* src, int bitsPerChannel)
{
    cerr << "This filter has not been implemented 7.\n";
    return NULL;
}

/* helper functions you may find useful for resampling */

/*
 * nearest neighbor sample
 */

Pixel ip_resample_nearest(Image* src, double x, double y)
{
    cerr << "This filter has not been implemented 8.\n";
    return Pixel(0,0,0);
}

/*
 * bilinear resample
 */

Pixel ip_resample_bilinear(Image* src, double x, double y)
{

    cerr << "This filter has not been implemented 9.\n";
    return Pixel(0,0,0);
}

/*
 * gausian sample
 */
Pixel ip_resample_gaussian(Image* src, double x, double y, int size, double sigma) 
{
    cerr << "This filter has not been implemented 10.\n";
    return Pixel(0,0,0);
}


/*
 * rotate image using one of three sampling techniques
 */
Image* ip_rotate (Image* src, double theta, int x, int y, int mode, 
                  int size, double sigma)
{
    cerr << "This filter has not been implemented 11.\n";
    return NULL;
}

/*
 * interpolate with the greyscale version of the image
 */
Image* ip_saturate (Image* src, double alpha)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* greyscale = ip_grey(src);
    Image* dest = new Image(width,height);
    dest = ip_interpolate (src, greyscale, alpha);
    return dest;
}


/*
 * scale image using one of three sampling techniques
 */
Image* ip_scale (Image* src, double xFac, double yFac, int mode, 
                 int size, double sigma)
{
    cerr << "This filter has not been implemented 12.\n";
    return NULL;
}

/*
 * create a new sepia tones image
 */
Image* ip_sepia (Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            double rp = 0.393*r + 0.769*g + 0.189*b;
            double gp = 0.349*r + 0.686*g + 0.168*b;
            double bp = 0.272*r + 0.534*g + 0.131*b;
            dest->setPixel_(w, h, 0, rp);
            dest->setPixel_(w, h, 1, gp);
            dest->setPixel_(w, h, 2, bp);
        }
    }
    return dest;
}


/*
 * create a new one bit/channel image with the intensity 0 or 1
 * depending on whether the input value is above or below the 
 * threshold
 */
Image* ip_threshold (Image* src, double cutoff)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* dest = new Image(width,height);
    for (int w=0; w<width; w++){
        for (int h=0;h<height; h++){
            double r = src->getPixel(w, h, 0);
            double g = src->getPixel(w, h, 1);
            double b = src->getPixel(w, h, 2);
            double rp = (r < cutoff) ? 0 : 1;
            double gp = (g < cutoff) ? 0 : 1;
            double bp = (b < cutoff) ? 0 : 1;
            dest->setPixel(w, h, 0, rp);
            dest->setPixel(w, h, 1, gp);
            dest->setPixel(w, h, 2, bp);
        }
    }
    return dest;
}




