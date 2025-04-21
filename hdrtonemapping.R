# Basic HDR tone mapping algorithm using contrast (sigmoid) curves
# www.overfitting.net
# https://www.overfitting.net/

library(tiff)  # save 16-bit TIFF's
library(terra)
library(png)  # save 8-bit PNG's


colourmatrix=function(img, colour=c(0.75, 0.5, 0.25), gamma=1) {
    # img must be a grayscale matrix
    
    # Gamma is applied before colouring
    img=replicate(3, img^(1/gamma))
    
    # Now the middle gray (0.5) becomes colour[]
    colour[colour<0.01]=0.01  # clip very low/high values
    colour[colour>0.99]=0.99
    colourgamma=log(0.5)/log(colour)
    for (i in 1:3) img[,,i]=img[,,i]^(1/colourgamma[i])

    return(img)
}

contrast=function(x, a=0.5, b=0.5, m=0, E=2) {
    f=log(b)/log(a)  # precalculate factor
    return(
        ifelse(  # vectorized ifelse
            x <= a,
            (m * x + (1 - m) * a * (x / a)^E) ^ f,
            (m * x + (1 - m) * (1 - (1 - a) * ((1 - x) / (1 - a))^E)) ^ f )
    )
}


#################################################
# 1. APPLY SIMPLE CONTRAST CURVE 

# Linear RAW development: dcraw -v -w -4 -T tiovivo.DNG
# resized to 1280x852 px
img=readTIFF("tiovivo.tif")^(1/2.2)  # read image and delinearize with 2.2 gamma
hist(img, breaks=800, xlim=c(0,1), ylim=c(0,14000))
a=median(img)  # turning point applied on median to maximize contrast split
b=a
m=0.1
E=2.5
abline(v=a, col='red')

imgcontrast=contrast(img, a, b, m, E)
hist(imgcontrast, breaks=800, xlim=c(0,1), ylim=c(0,14000))
acontrast=median(imgcontrast)
paste0("Before/after median is preserved: ",
       round(a,4), " -> ", round(acontrast,4))
abline(v=acontrast, col='red')
writeTIFF(imgcontrast, "tiovivocontrast.tif", bits.per.sample=16)






