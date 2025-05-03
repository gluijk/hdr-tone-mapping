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
    if (a==0 | a==1) {  # if x=0 or x=1 the identity function is applied
        return(x)
    } else {
        f=log(b)/log(a)  # precalculate factor
        return(
            ifelse(  # vectorized ifelse
                x <= a,
                (m * x + (1 - m) * a * (x / a)^E) ^ f,
                (m * x + (1 - m) * (1 - (1 - a) * ((1 - x) / (1 - a))^E)) ^ f )
        )
    }
}

# Resample
arrayresample=function(img, DIMX, DIMY, method='bilinear') {
    # method=c('near', 'bilinear', 'cubic', 'cubicspline', 'lanczos')
    
    require(terra)
    
    raster=rast(img)
    rasterrs=rast(nrows=round(DIMY), ncols=round(DIMX), extent=ext(raster))
    rasterrs=resample(raster, rasterrs, method=method, threads=TRUE)
    
    if (is.matrix(img)) return (matrix(as.array(rasterrs), nrow=nrow(rasterrs)))    
    return (as.array(rasterrs))  # convert back to matrix/array
}


# Blur
# https://stackoverflow.com/questions/70429190/how-can-i-perform-neighborhood-analysis-in-terra-or-raster-and-keep-the-same-na
arrayblur=function(img, radius=10, kernel='gaussian', fun='mean') {
    # radius: radius of the averaging window
    # kernel: 'gaussian', 'spherical'. Otherwise all 1's
    
    require(terra)
    
    # Build circular kernel
    D=2*radius+1  # D will always be an odd number as required by focal()
    sigma=radius/3  # keep gaussian shape for every radius
    w=matrix(1, nrow=D, ncol=D)
    if (kernel=='gaussian') {  # gaussian filter
        w=exp(-((row(w)-(radius+1))^2 + (col(w)-(radius+1))^2) / (2 * sigma^2))
    } else if (kernel=='spherical') {  # spherical filter
        w=1 - ((row(w)-(radius+1))^2 + (col(w)-(radius+1))^2) / (radius+1)^2        
    }
    # Ignore values out of the radius
    w[(row(w)-(radius+1))^2 + (col(w)-(radius+1))^2 > (radius+1)^2]=NA
    writePNG(w, "blurkernel.png")
    
    raster=rast(img)  # process as raster
    rasterblur=focal(raster, w=w, fun=fun, na.rm=TRUE, na.policy='omit')
    
    if (is.matrix(img)) return (matrix(as.array(rasterblur), nrow=nrow(rasterblur)))
    else return (as.array(rasterblur))  # convert back to matrix/array
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


#################################################
# 2. BASIC TONE MAPPING

# Create blurred version
DEM=readTIFF("tiovivo.tif")^(1/2.2)  # read image and delinearize with 2.2 gamma
DEM=readTIFF("room.tif")  # read image
DEM=readTIFF("peninsuladem.tif")  # read image
# DEMblur=0.299*DEM[,,1]+0.587*DEM[,,2]+0.114*DEM[,,3]  # B&W blur
DEMblur=DEM
DEMblur=arrayblur(DEMblur, radius=80)
writeTIFF(DEMblur, "DEMblur.tif", bits.per.sample=16)

# Tone mapping (looped version)

# Reduce global contrast (static curve 1)
a1=0.5  # median(DEM[DEM>0 & DEM<1])  # 0.5
b1=0.5
m1=0.2
E1=0.2  # 0.7
# Increase local contrast (adaptive curve 2)
m2=0.3
E2=4  # 1.2
DEMtonemap=DEM*0


# Draw Tone mapping curve set
png("tonemappingcurveset.png", width=800, height=800)
x=seq(0, 1, 0.001)
y=contrast(x, a1, b1, m1, E1)
plot(x, y, type='l', asp=1, xlim=c(0,1), ylim=c(0,1),
     xaxt="n", yaxt="n",
     main=paste0('Tone mapping curve set\n(a1=', a1, ', b1=', b1, ', m1=', m1,
                 ', E1=', E1, ' / m2=', m2, ', E2=', E2, ')'),
            xlab="IN", ylab="OUT")
lines(c(0,1), c(0,1), type='l', col='gray')
N=21
for (i in seq(0, 1, length.out=N)[2:(N-1)]) {
    a2=i
    b2=contrast(a2, a1, b1, m1, E1)
    y=contrast(x, a2, b2, m2, E2)
    lines(x, y, type='l', col='red', asp=1,
          xlim=c(0,1), ylim=c(0,1))
    points(a2, b2, pch=16, col="red")
}
abline(h=c(0,0.5,1), v=c(0,0.5,1), lty=2, col='gray')
axis(1, at=c(0, 0.5, 1))
axis(2, at=c(0, 0.5, 1))
dev.off()


# If colour image
for (i in 1:nrow(DEM)) {
    for (j in 1:ncol(DEM)) {
        a2=DEMblur[i,j]
        b2=contrast(a2, a1, b1, m1, E1)
        DEMtonemap[i,j,]=contrast(DEM[i,j,], a2, b2, m2, E2)
    }
}
writeTIFF(DEMtonemap, "DEMtonemap_room.tif", bits.per.sample=16)


# If B&W image
for (i in 1:nrow(DEM)) {
    for (j in 1:ncol(DEM)) {
        a2=DEMblur[i,j]
        b2=contrast(a2, a1, b1, m1, E1)
        DEMtonemap[i,j]=contrast(DEM[i,j], a2, b2, m2, E2)
    }
}
writeTIFF(DEMtonemap, "DEMtonemap_DEM.tif", bits.per.sample=16)
