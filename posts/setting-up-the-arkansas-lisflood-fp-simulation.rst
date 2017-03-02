.. title: Setting up the LISFLOOD-FP Arkansas simulation
.. slug: Setting-up-the-LISFLOOD-FP-Arkansas-simulation
.. date: 2016-12-27 13:43:00
.. tags: lisflood, swot
.. link: 
.. description: 
.. type: text
.. author: Kostas Andreadis

1 Setting up the Arkansas LISFLOOD-FP simulation
------------------------------------------------

We start by downloading the HydroSHEDS `elevation DEM <http://earlywarning.usgs.gov/hydrodata/sa_30s_zip_grid/na_dem_30s_grid.zip>`_ and reprojecting it to the World Mercator projection. Additionally we download shapefiles for the basin and rivers from `NOHRSC <http://www.nohrsc.noaa.gov/gisdatasets/>`_.

.. TEASER_END

.. code:: shell

    cd ~/Work/Swot
    gdalwarp -dstnodata -9999 -t_srs epsg:3395 data/hydrosheds/na_dem_15s/na_dem_15s na_dem.tif

Then we extract the extent of the basin shapefile and use those coordinates as the bounding box to generate the Arkansas River basin DEM.

.. code:: shell

    cd ~/Work/Swot
    ogrinfo -al -so ~/Downloads/b_abrfc/arkansas.shp | grep Extent | awk -F'[(,]' '{printf("%f %f\n%f %f\n",$2,$5,$4,$3)}' | gdaltransform -s_srs epsg:4326 -t_srs epsg:3395 -output_xy
    gdal_translate -a_nodata -9999 -projwin -11866240.2704725 4749554.46189774 -10240048.0795739 3922731.2771146 na_dem.tif data/hydrosheds/arkansas_dem.tif

Within GRASS GIS, we calculate the flow accumulation and derive a consistent river network using 10,000 km\ :sup:`2`\ as a drainage area threshold. Before deriving the river network we create a depressionless DEM (i.e. hydrologic conditioning).

.. code:: shell

    cd ~/Work/Swot
    v.in.ogr in=data/arkansas_basin.shp out=basin
    v.to.rast in=basin out=basin use=val value=1
    r.in.gdal in=data/hydrosheds/arkansas_dem.tif out=hydrosheds.dem
    g.region rast=hydrosheds.dem
    g.region res=1000
    r.resamp.stats in=hydrosheds.dem out=hydrosheds.elev method=minimum
    r.hydrodem -a in=hydrosheds.elev out=hydrosheds.felev
    r.watershed -s elev=hydrosheds.felev accum=hydrosheds.flowacc stream=hydrosheds.river drain=hydrosheds.flowdir threshold=10000
    r.out.gdal in=hydrosheds.flowacc out=data/hydrosheds/arkansas_acc.tif nodata=-9999
    r.out.gdal in=hydrosheds.flowdir out=data/hydrosheds/arkansas_flowdir.tif nodata=-9999
    r.out.gdal type=Float64 in=hydrosheds.felev out=data/hydrosheds/arkansas_elev.tif nodata=-9999

Then we import the width and depth `database <https://zenodo.org/record/61758#.WF8A57YrKRs>`_ and attach the attributes to the nearest river segment in the derived network.

.. code:: shell

    cd ~/Work/Swot
    ogr2ogr -t_srs epsg:3395 -clipdst -11866240.2704725 3922731.2771146 -10240048.0795739 4749554.46189774 rivs.shp data/hydrosheds/nariv.shp
    v.in.ogr -r in=rivs.shp out=nariv where='AREA>10000'
    r.mask basin
    r.to.vect in=hydrosheds.river out=river type=line -s
    r.mask -r
    v.db.addcolumn map=river col='width real'
    v.db.addcolumn map=river col='depth real'
    v.distance from=river to=nariv upload=to_attr to_column=DEPTH column=depth
    v.distance from=river to=nariv upload=to_attr to_column=WIDTH column=width

Now we rasterize the river vector to get the width and depth rasters.

.. code:: shell

    cd ~/Work/Swot
    v.to.rast in=river out=hydrosheds.width use=attr attribute_column=width
    v.to.rast in=river out=hydrosheds.depth use=attr attribute_column=depth
    r.out.gdal in=hydrosheds.width out=data/hydrosheds/arkansas_widths.tif nodata=-9999
    r.out.gdal in=hydrosheds.depth out=data/hydrosheds/arkansas_depths.tif nodata=-9999

The next steps involve preparing the actual LISFLOOD-FP input files. Given the uncertainty in the DEM, we will smooth the river bank heights to avoid any numerical instabilities during simulation. 
This is accomplished by first identifying the upstream and downstream boundary points in the domain, labeling the river channels and calculating each one's chainage (i.e. distance downstream). 

.. code:: ipython

    import sys
    sys.path.append("/Users/kandread/Work/Swot/scripts")
    import lisflood
    flowdir = lisflood.read_raster("/Users/kandread/Work/Swot/data/hydrosheds/arkansas_flowdir.tif")
    flowacc = lisflood.read_raster("/Users/kandread/Work/Swot/data/hydrosheds/arkansas_acc.tif")
    widths = lisflood.read_raster("/Users/kandread/Work/Swot/data/hydrosheds/arkansas_widths.tif")
    river = widths > 0
    bndpts = lisflood.find_boundary_points(river, flowdir)
    labels, chainage = lisflood.calc_chainage(river, flowdir, abs(flowacc), bndpts, 1000)

The smoothing of the bank heights is done by using a `LOWESS <http://statsmodels.sourceforge.net/devel/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html>`_ local regression, and the channel is burned in to the DEM by subtracting the depth raster.

.. code:: ipython

    elev = lisflood.read_raster("/Users/kandread/Work/Swot/data/hydrosheds/arkansas_elev.tif")
    selev = lisflood.smooth_bank_heights(elev, labels, chainage)
    lisflood.write_raster(selev, "/Users/kandread/Work/Swot/data/hydrosheds/arkansas_selev.tif", "/Users/kandread/Work/Swot/data/hydrosheds/arkansas_elev.tif")
    depths = lisflood.read_raster("/Users/kandread/Work/Swot/data/hydrosheds/arkansas_depths.tif")
    depths[depths < 0] = 0.0
    belev = selev - depths
    lisflood.write_raster(belev, "/Users/kandread/Work/Swot/data/hydrosheds/arkansas_bed.tif", "/Users/kandread/Work/Swot/data/hydrosheds/arkansas_elev.tif")
    y0 = 0.8 * depths
    y0[widths <= 0] = 0.0
    lisflood.write_raster(y0, "/Users/kandread/Work/Swot/data/hydrosheds/arkansas_initwd.tif", "/Users/kandread/Work/Swot/data/hydrosheds/arkansas_elev.tif")

Next we identify the upstream and lateral inflow points, and generate the BCI file.

.. code:: ipython

    from osgeo import gdal
    inflows = lisflood.identify_inflows(river, chainage, labels, abs(flowacc), 10000)
    f = gdal.Open("/Users/kandread/Work/Swot/data/hydrosheds/arkansas_selev.tif")
    xul, xres, _, yul, _, yres = f.GetGeoTransform()
    f = None
    nrows, ncols = elev.shape
    lisflood.write_bci("/Users/kandread/Work/Swot/input/arkansas.bci", inflows, xul, yul, xres, yres, nrows, ncols)

Then we generate the DEM and sub-grid channel width files.

.. code:: shell

    cd ~/Work/Swot
    gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 data/hydrosheds/arkansas_selev.tif input/arkansas.dem
    gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 data/hydrosheds/arkansas_widths.tif input/arkansas.width
    gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 data/hydrosheds/arkansas_bed.tif input/arkansas.bed
    gdal_translate -tr 1000 -1000 -of AAIGrid -a_nodata -9999 -co DECIMAL_PRECISION=2 data/hydrosheds/arkansas_initwd.tif input/arkansas.initwd

If we need LISFLOOD-FP to produce a time series of river flow at specific locations, we need to generate a virtual gauge file.

.. code:: ipython

    x = [-10646496, -10281401]
    y = [4251469, 4111473]
    lisflood.write_gauges(x, y, flowdir, "/Users/kandread/Work/Swot/data/hydrosheds/arkansas_widths.tif", "/Users/kandread/Work/Swot/input/arkansas.gauge")

The final step is generating the BDY file for LISFLOOD-FP, which contains the streamflow at the inflow points (in m\ :sup:`2`\/s). For this simulation we will use the output of the `VIC routing model <https://github.com/UW-Hydro/VIC_Routing>`_ forced by the `NLDAS-2 VIC model output <http://ldas.gsfc.nasa.gov/nldas/NLDAS2model.php>`_. We need to prepare the input files for the routing model, describing the flow direction, flow fraction and station locations. 
The station file is created from the inflows identified

.. code:: ipython

    from pyproj import Proj
    wmerc = Proj("+init=EPSG:3395")
    fout = open("arkansas.stations", 'w')
    for sta, xy in enumerate(inflows):
        x, y = wmerc(xul+xres*xy[1], yul+yres*xy[0], inverse=True)
        xi = int((x+106.625)/0.125) + 1
        yi = int((y-33.375)/0.125) + 1
        fout.write("1\tAR{0:03d}\t{1}\t{2}\t-9999\nNONE\n".format(sta, xi, yi))
    fout.close()

Within GRASS GIS, we create a new region with Lat/Long projection and work within that region to create the necessary files.

.. code:: shell

    g.proj epsg=4326 location=arkansas.latlon
    g.mapset mapset=PERMANENT location=arkansas.latlon
    g.region `r.proj location=arkansas mapset=arkansas in=hydrosheds.elev -g`
    #g.region n=39.5 s=32.25 e=-91.875 w=-106.75 res=0:00:30
    r.proj mapset=arkansas location=arkansas in=hydrosheds.felev method=bilinear
    r.proj mapset=arkansas location=arkansas in=basin
    r.mapcalc exp='basin1=if(isnull(basin),0,1)'
    g.region res=0.125 -a
    r.watershed -s elev=hydrosheds.felev accum=hydrosheds.flowacc drain=hydrosheds.flowdir stream=hydrosheds.river threshold=15
    r.resamp.stats in=basin1 out=fract
    r.reclass in=hydrosheds.flowdir out=flowdir rules=flowdir.rules
    r.mapcalc --o exp='flowdir=flowdir'
    r.null map=flowdir setnull=0
    r.null map=fract null=0
    r.out.gdal in=fract out=arkansas.fract format=AAIGrid nodata=0
    r.out.gdal in=flowdir out=arkansas.flowdir format=AAIGrid nodata=0
    awk 'BEGIN{OFS="|"}/^1/{print(-106.625+0.125*$3,33.375+0.125*$4,$2)}' arkansas.stations | v.in.ascii in=- out=inflows col='x real,y real,name text'
    r.stream.snap in=inflows out=stations accum=hydrosheds.flowacc stream=hydrosheds.river radius=3
    v.to.rast in=stations out=inflows use=cat
    r.stats in=inflows -n -x | sort -k3 -n | awk '{printf("1\tAR%03d\t%d\t%d\t-9999\nNONE\n",$3,$1,49-$2+1)}' > arkansas.stations
    v.db.addtable map=stations
    v.db.addcolumn map=stations col='flowacc real'
    v.what.rast map=stations raster=hydrosheds.flowacc column=flowacc
    v.db.select stations sep=" " | sort -k1 -n | awk 'function abs(x){return ((x < 0.0) ? -x : x)}NR>1{print(abs($2)*157.828125)}' > stations.flowacc

Finally, the BDY file is written using the VIC routing model's output.

.. code:: ipython

    import numpy as np
    stations = ["AR{0:03d}".format(s+1) for s in range(len(inflows))]
    a0 = np.array([abs(flowacc[i[0],i[1]]) for i in inflows])
    a = np.loadtxt("stations.flowacc")
    area_mult = a0 / a
    lisflood.write_bdy("/Users/kandread/Work/Swot/input/arkansas.bdy", "/Volumes/External2/nldas2", stations, area_mult)
