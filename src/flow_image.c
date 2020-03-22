#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"

// Draws a line on an image with color corresponding to the direction of line
// image im: image to draw line on
// float x, y: starting point of line
// float dx, dy: vector corresponding to line angle and magnitude
void draw_line(image im, float x, float y, float dx, float dy)
{
    assert(im.c == 3);
    float angle = 6*(atan2(dy, dx) / TWOPI + .5);
    int index = floor(angle);
    float f = angle - index;
    float r, g, b;
    if(index == 0){
        r = 1; g = f; b = 0;
    } else if(index == 1){
        r = 1-f; g = 1; b = 0;
    } else if(index == 2){
        r = 0; g = 1; b = f;
    } else if(index == 3){
        r = 0; g = 1-f; b = 1;
    } else if(index == 4){
        r = f; g = 0; b = 1;
    } else {
        r = 1; g = 0; b = 1-f;
    }
    float i;
    float d = sqrt(dx*dx + dy*dy);
    for(i = 0; i < d; i += 1){
        int xi = x + dx*i/d;
        int yi = y + dy*i/d;
        set_pixel(im, xi, yi, 0, r);
        set_pixel(im, xi, yi, 1, g);
        set_pixel(im, xi, yi, 2, b);
    }
}

// Make an integral image or summed area table from an image
// image im: image to process
// returns: image I such that I[x,y] = sum{i<=x, j<=y}(im[i,j])
image make_integral_image(image im)
{
    image integ = make_image(im.w, im.h, im.c);
    float im_pix; //original image pixel @i(x,y) 
    float integ_y1 = 0.0; //integ imag @I(x,y-1)
    float integ_x1 = 0.0; //integ imag @I(x-1,y)
    float integ_x1_y1 = 0.0; //integ imag @I(x-1,y-1)
    float new_pix = 0.0; //new integ pixel @I(x,y)

    for(int c = 0; c < im.c; ++c)
    {
        for(int x = 0; x < im.w; ++x)
        {   
          
            for(int y = 0; y < im.h; ++y)
            {
                im_pix = get_pixel(im,x,y,c);
                integ_y1 = 0.0;
                integ_x1 = 0.0;
                integ_x1_y1 = 0.0;

                if(y>0)
                {
                    integ_y1 = get_pixel(integ,x,y-1,c); 
                }
                
                
                if(x>0)
                {
                    integ_x1 = get_pixel(integ,x-1,y,c); 
                }
                
                if(y>0 && x >0)
                {
                    integ_x1_y1 = get_pixel(integ,x-1,y-1,c); 
                }
                

                new_pix = im_pix + integ_y1 + integ_x1 - integ_x1_y1; 
                set_pixel(integ,x,y,c,new_pix); 
            }

        }
    }
    return integ;
}

// Apply a box filter to an image using an integral image for speed
// image im: image to smooth
// int s: window size for box filter
// returns: smoothed image

image box_filter_image(image im, int s)
{
    int i,j,k;
    image integ = make_integral_image(im);
    image S = make_image(im.w, im.h, im.c);
    // TODO: fill in S using the integral image.
    int t; //top pixel 
    int r; //right pixel 
    int l; //left pixel 
    int b; //bottom pixel 
    int filter_w = s/2; 
    float b_l; //bottom left
    float b_r; //bottom right
    float t_l; //top left
    float t_r; //top right


    for(int c = 0; c < integ.c; ++c)
    {
        for(int x = 0; x < integ.w; ++x)
        {   
            l = x - filter_w -1; 
            r = x + filter_w; 
            for(int y = 0; y < integ.h; ++y)
            {   
                t = y - filter_w -1; 
                b = y + filter_w; 

                if(t < 0 && l >= 0)
                {
                    t_r = 0.; 
                    t_l = 0.; 
                    b_l = get_pixel(integ,l,b,c); 
                }
                else if (l < 0 && t >= 0)
                {   
                    b_l = 0.; 
                    t_l = 0.; 
                    t_r = get_pixel(integ,r,t,c); 

                }
                else if(l < 0 && t < 0)
                {
                    t_l = 0.; 
                    t_r = 0.;
                    b_l = 0.; 
                }
                else 
                {
                    t_l = get_pixel(integ,l,t,c); 
                    t_r = get_pixel(integ,r,t,c); 
                    b_l = get_pixel(integ,l,b,c);
                }
                
                b_r = get_pixel(integ,r,b,c); 
                float sum = b_r + t_l - b_l - t_r; 
                set_pixel(S,x,y,c,sum); 
            }

        }
    }
    free_image(integ); 
    return S;
}
// Calculate the time-structure matrix of an image pair.
// image im: the input image.
// image prev: the previous image in sequence.
// int s: window size for smoothing.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          3rd channel is IxIy, 4th channel is IxIt, 5th channel is IyIt.
image time_structure_matrix(image im, image prev, int s)
{
    int i;
    int converted = 0;
    image S = make_image(im.w,im.h,5); 
    if(im.c == 3){
        converted = 1;
        im = rgb_to_grayscale(im);
        prev = rgb_to_grayscale(prev);
    }

    // TODO: calculate gradients, structure components, and smooth them
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();

    image I_x = convolve_image(im, gx_filter, 0);
    image I_y = convolve_image(im, gy_filter, 0);

    for(int x = 0; x < S.w; ++x)
    {
        for(int y = 0; y < S.h; ++y)
        {
            float pix_Ix = get_pixel(I_x,x,y,0);  
            float pix_Iy = get_pixel(I_y,x,y,0);
            float pix_It = get_pixel(im,x,y,0) - get_pixel(prev,x,y,0);   //only using grayscale 
            set_pixel(S,x,y,0,pix_Ix * pix_Ix);
            set_pixel(S,x,y,1,pix_Iy * pix_Iy);
            set_pixel(S,x,y,2,pix_Ix * pix_Iy);
            set_pixel(S,x,y,3,pix_Ix * pix_It);
            set_pixel(S,x,y,4,pix_Iy * pix_It);

        }
    }

    S = box_filter_image(S,s); 
    free_image(I_x); 
    free_image(I_y);
    free_image(gx_filter); 
    free_image(gy_filter); 

    if(converted){
        free_image(im); free_image(prev);
    }
    return S;
}
// Calculate the velocity given a structure image
// image S: time-structure image
// int stride: only calculate subset of pixels for speed
image velocity_image(image S, int stride)
{
    image v = make_image(S.w/stride, S.h/stride, 3);
    int i, j;
    matrix M = make_matrix(2,2);
    matrix Mt = make_matrix(2,1); 
    for(j = (stride-1)/2; j < S.h; j += stride){
        for(i = (stride-1)/2; i < S.w; i += stride){
            float Ixx = S.data[i + S.w*j + 0*S.w*S.h];
            float Iyy = S.data[i + S.w*j + 1*S.w*S.h];
            float Ixy = S.data[i + S.w*j + 2*S.w*S.h];
            float Ixt = S.data[i + S.w*j + 3*S.w*S.h];
            float Iyt = S.data[i + S.w*j + 4*S.w*S.h];

            // TODO: calculate vx and vy using the flow equation
            float vx = 0;
            float vy = 0;
            M.data[0][0] = Ixx;
            M.data[1][0] = Ixy; 
            M.data[0][1] = Ixy; 
            M.data[1][1] = Iyy;

            Mt.data[0][0] = -Ixt; 
            Mt.data[1][0] = -Iyt; 

            matrix M_inverse = matrix_invert(M);

            if(M_inverse.rows > 0)
            {
                matrix Vel = matrix_mult_matrix(M_inverse,Mt);
                vx = Vel.data[0][0]; 
                vy = Vel.data[1][0];  
                free_matrix(Vel); 
            }
            set_pixel(v, i/stride, j/stride, 0, vx);
            set_pixel(v, i/stride, j/stride, 1, vy);
            free_matrix(M_inverse); 
        }
    }
    free_matrix(M);
    free_matrix(Mt);
    return v;
}

// Draw lines on an image given the velocity
// image im: image to draw on
// image v: velocity of each pixel
// float scale: scalar to multiply velocity by for drawing
void draw_flow(image im, image v, float scale)
{
    int stride = im.w / v.w;
    int i,j;
    for (j = (stride-1)/2; j < im.h; j += stride) {
        for (i = (stride-1)/2; i < im.w; i += stride) {
            float dx = scale*get_pixel(v, i/stride, j/stride, 0);
            float dy = scale*get_pixel(v, i/stride, j/stride, 1);
            if(fabs(dx) > im.w) dx = 0;
            if(fabs(dy) > im.h) dy = 0;
            draw_line(im, i, j, dx, dy);
        }
    }
}


// Constrain the absolute value of each image pixel
// image im: image to constrain
// float v: each pixel will be in range [-v, v]
void constrain_image(image im, float v)
{
    int i;
    for(i = 0; i < im.w*im.h*im.c; ++i){
        if (im.data[i] < -v) im.data[i] = -v;
        if (im.data[i] >  v) im.data[i] =  v;
    }
}

// Calculate the optical flow between two images
// image im: current image
// image prev: previous image
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// returns: velocity matrix
image optical_flow_images(image im, image prev, int smooth, int stride)
{
    image S = time_structure_matrix(im, prev, smooth);   
    image v = velocity_image(S, stride);
    constrain_image(v, 6);
    image vs = smooth_image(v, 2);
    free_image(v);
    free_image(S);
    return vs;
}

// Run optical flow demo on webcam
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// int div: downsampling factor for images from webcam
void optical_flow_webcam(int smooth, int stride, int div)
{
#ifdef OPENCV
    CvCapture * cap;
    cap = cvCaptureFromCAM(0);
    image prev = get_image_from_stream(cap);
    image prev_c = nn_resize(prev, prev.w/div, prev.h/div);
    image im = get_image_from_stream(cap);
    image im_c = nn_resize(im, im.w/div, im.h/div);
    while(im.data){
        image copy = copy_image(im);
        image v = optical_flow_images(im_c, prev_c, smooth, stride);
        draw_flow(copy, v, smooth*div);
        int key = show_image(copy, "flow", 5);
        free_image(v);
        free_image(copy);
        free_image(prev);
        free_image(prev_c);
        prev = im;
        prev_c = im_c;
        if(key != -1) {
            key = key % 256;
            printf("%d\n", key);
            if (key == 27) break;
        }
        im = get_image_from_stream(cap);
        im_c = nn_resize(im, im.w/div, im.h/div);
    }
#else
    fprintf(stderr, "Must compile with OpenCV\n");
#endif
}