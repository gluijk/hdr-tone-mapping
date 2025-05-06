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


solid=function(DEM, altitude=0, isnan=0) {
    # solid() calculates a solid version of a DEM
    #
    # altitude: DEM altitude level contour
    # isnan: value assigned to NaN data
    
    DIMY=nrow(DEM)
    DIMX=ncol(DEM)
    
    # If array turn its first dimension into a genuine matrix
    if (!is.matrix(DEM)) {
        print("WARNING: input DEM is not a matrix but an array. First dimension is used")
        DEM=matrix(DEM[,,1], nrow=DIMY, ncol=DIMX)
    }
    DEM[is.nan(DEM)]=isnan
    
    # Calculate solid map from DEM
    solidmap=DEM*0
    solidmap[DEM > altitude]=1
    
    return(solidmap)
}


#################################################
# 1. SIMPLE CONTRAST CURVE 

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



###########################################################
# ENHANCED DEM (REM) THROUGH HDR TONE MAPPING

###########################################################

# 2. READ DEM 

# https://data.humdata.org/dataset/worldpop-population-density-for-spain
peninsula=rast("gebco_2024_n45.8569_s34.2114_w-11.5137_e5.2295.tif")
peninsula
plot(peninsula)
RESOLUTION=res(peninsula)[1]  # 0.004167 degrees resolution

# REPROJECT raster from Longitude Latitude (+proj=longlat)/WGS84
# to Lambert Conic Conformal (+proj=lcc)/WGS84
# by default crs="+proj=lcc +ellps=WGS84 +lat_1=33 +lat_2=45 +lon_0=0"
# CRS="+proj=lcc +ellps=WGS84 +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=0 +units=km"  # default
# CRS="+proj=lcc +ellps=WGS84 +lat_1=37 +lat_2=43 +lat_0=40 +lon_0=-3 +units=km"  # better for Peninsula
# CRS="+proj=lcc +ellps=WGS84 +lat_1=30 +lat_2=44 +lat_0=37 +lon_0=-3 +units=km"  # better for Spain incl. Canarias
# CRS="+proj=merc +ellps=WGS84 +datum=WGS84 +units=km"  # Mercator for Spain incl. Canarias
CRS="+proj=lcc +ellps=WGS84 +lat_1=37 +lat_2=43 +lat_0=40 +lon_0=-3.141667 +units=km"  # better for Peninsula

peninsula=project(x=peninsula, y=CRS, threads=TRUE)
peninsula
plot(peninsula)
abline(h=0, v=0)  # centre of reference
RESOLUTION=res(peninsula)[1]  # 0.392 km resolution

# CROP raster to drop NaN's
cropdef=ext(-659, 659, -600, 600)
peninsulacrop=crop(x=peninsula, y=cropdef)
peninsulacrop
plot(peninsulacrop)
abline(h=0, v=0)  # centre of reference

# RESAMPLE raster to Full HD
# DIMY=1080
DIMX=1920
DIMY=round(DIMX*nrow(peninsulacrop)/ncol(peninsulacrop))  # 1748 px
peninsulars=rast(nrows=DIMY, ncols=DIMX, extent=ext(peninsulacrop))
peninsulars=resample(x=peninsulacrop, y=peninsulars, method='bilinear', threads=TRUE)
plot(peninsulars)
abline(h=0, v=0)  # centre of reference
RESOLUTION=res(peninsulars)[1]  # 0.686 km resolution


###########################################################

# 3. CONVERT TO MATRIX AND CALCULATE SOLID AND CONTOUR

# Convert to matrix and save as TIFF
# DEM=matrix(as.array(peninsulars), nrow=nrow(peninsulars))
DEM=as.matrix(peninsulars, wide=TRUE)  # no need to use as.array()
hist(DEM, breaks=1000)
DEM[DEM < 0]=0
writeTIFF((DEM/max(DEM))^(1/2.2), "peninsuladem.tif", bits.per.sample=16, compression='LZW')

# Solid map for masking in Photoshop
DEMsolid=solid(DEM)
writeTIFF(DEMsolid, "peninsulasolid.tif", compression='LZW')


#################################################
# 4. BASIC TONE MAPPING

# Create blurred version
name="peninsuladem"  # "room", "tiovivo"
gamma=1  # 2.2 for RAW development
DEM=readTIFF(paste0(name, ".tif"))^(1/gamma)  # read image
if (length(dim(DEM))==2) {  # B&W image
        DEMblur=DEM
    } else {  # colour image
        DEMblur=0.299*DEM[,,1]+0.587*DEM[,,2]+0.114*DEM[,,3]  # B&W blur
    }
DEMblur=arrayblur(DEMblur, radius=80)  # takes time...
writeTIFF(DEMblur, paste0("DEMblur_", name, ".tif"), bits.per.sample=16)


# Tone mapping (looped version)

# Reduce global contrast (static curve 1)
a1=0.5  # median(DEM[DEM>0 & DEM<1])  # 0.5
b1=0.5
m1=0
E1=0  # 0.7
# Increase local contrast (adaptive curve 2)
m2=0
E2=5  # 1.2


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
writeTIFF(DEMtonemap, paste0("DEMtonemap_", name, ".tif"), bits.per.sample=16)


#################################################
# 5. 3D MESHES

z <- readTIFF("peninsuladem.tif")^2.2  # undo gamma 2.2
z <- readTIFF("DEMtonemap_peninsuladem.tif")^2.2  # undo gamma 2.2
z[z==0]=NA

# Create x and y coordinate sequences
x <- 1:nrow(z)
y <- 1:ncol(z)

# Create 3D mesh plot
persp(x, y, z,
      theta = 90, phi = 45,      # viewing angles
      expand = 0.01,  # 0.5,              # vertical exaggeration
      col = "lightblue",         # color of surface
      shade = 0.5,               # shading
      border = NA,               # remove borders for a smooth look
      axes = FALSE,
      xlab = "", ylab = "", zlab = "")     # axes ticks

