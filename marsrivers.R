# Basic HDR tone mapping algorithm using contrast (sigmoid) curves
# www.overfitting.net
# https://www.overfitting.net/

library(tiff)  # save 16-bit TIFF's
library(terra)
library(png)  # save 8-bit PNG's


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


###########################################################
# ENHANCED DEM (REM) THROUGH HDR TONE MAPPING


# 1. READ DEM 

# READ raster
mars=rast("Mars_MGS_MOLA_DEM_mosaic_global_463m.tif")  # read GeoTIFF file
mars
# dimensions  : 23040, 46080 pixels
# resolution  : 463.0935 m/pixel
# extent      : -10669675, 10669675, -5334838, 5334838 m
# coord. ref. : Equirectangular Mars
plot(mars)
RESOLUTION=res(mars)[1]


# CROP raster to area of interest
CENTREX=-3.2e6
CENTREY=1e6
DELTADISTX=8e6*0.65
# DELTADISTY=round(DELTADISTX*1080/1920)
DELTADISTY=round(DELTADISTX*1920/1920)
cropdef=ext(CENTREX-DELTADISTX/2, CENTREX+DELTADISTX/2,
            CENTREY-DELTADISTY/2, CENTREY+DELTADISTY/2)
marscrop=crop(x=mars, y=cropdef)
marscrop
plot(marscrop)


# RESAMPLE raster to Full HD
DIMY=1920
DIMX=1920
marscroprs=rast(nrows=DIMY, ncols=DIMX, extent=ext(marscrop))
marscroprs=resample(x=marscrop, y=marscroprs, method='bilinear', threads=TRUE)
marscroprs
plot(marscroprs)
RESOLUTION=res(marscroprs)[1]


###########################################################

# 2. CONVERT TO MATRIX AND CALCULATE SOLID AND CONTOUR

# Convert to matrix and save as TIFF
DEM=as.matrix(marscroprs, wide=TRUE)
hist(DEM, breaks=800)
DEM=DEM-min(DEM)
DEM=DEM/max(DEM)
writeTIFF(DEM, "mars.tif", bits.per.sample=16, compression='LZW')


#################################################
# 3. BASIC TONE MAPPING

# Create blurred version
if (length(dim(DEM))==2) {  # B&W image
        DEMblur=DEM
    } else {  # colour image
        DEMblur=0.299*DEM[,,1]+0.587*DEM[,,2]+0.114*DEM[,,3]  # B&W blur
    }
DEMblur=arrayblur(DEMblur, radius=80)  # takes time...
writeTIFF(DEMblur, "DEMblurmars.tif", bits.per.sample=16)


# Tone mapping (looped version)

# Reduce global contrast (static curve 1)
a1=0.5  # median(DEM[DEM>0 & DEM<1])  # 0.5
b1=0.5
m1=0
E1=0.3  # 0  # 0.7
# Increase local contrast (adaptive curve 2)
m2=0
E2=3  # 1.2


# Draw Tone mapping curve set
png("tonemappingcurveset.png", width=512, height=512)
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


# Generate output tone mapped image
DEMtonemap=DEM*0
if (length(dim(DEM))==2) {  # B&W image
    for (i in 1:nrow(DEM)) {
        for (j in 1:ncol(DEM)) {
            a2=DEMblur[i,j]
            b2=contrast(a2, a1, b1, m1, E1)
            DEMtonemap[i,j]=contrast(DEM[i,j], a2, b2, m2, E2)
        }
    }
} else {  # colour image
    for (i in 1:nrow(DEM)) {
        for (j in 1:ncol(DEM)) {
            a2=DEMblur[i,j]
            b2=contrast(a2, a1, b1, m1, E1)
            DEMtonemap[i,j,]=contrast(DEM[i,j,], a2, b2, m2, E2)
        }
    }
}
writeTIFF(DEMtonemap, "DEMtonemapmars.tif", bits.per.sample=16)
