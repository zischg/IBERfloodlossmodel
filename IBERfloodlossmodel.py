#******************************************************************************
#Script for attributing flow deoths from IBER simulations to building polygons
#Andreas Paul Zischg, 12.08.2017
#Requires:
#IBER mesh file .msh (triangles only)
#IBER results file .res (with flow depths on nodes per timestep)
#buildings polygon shapefile with sepcified columns:
#   Wert ... Reconstruction value of building
#   AnzPers ... number of residents
#   Orig_FID ... unique feature identifier
#******************************************************************************

import arcpy
import numpy
import math
import os
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

arcpy.env.workspace = myworkspace
arcpy.env.overwriteOutput = True

#set environment and workspace
myworkspace = "/home"
tempdir="/temp"
outputdelimiter="\t"
flowdepthstatmethod="MEAN"#"MAX" or "MEAN"
searchradius=0.5
gebshp="\home\buildings.shp" #shapefile of buildings with the obligatory column names:
#"Orig_FID" .. a unique ID (integer) of each building, 
#"Wert" .. "the monetary value of the building (DOUBLE), 
#"AnzPers" .. the number of residents in the building (integer) 
floodplainname="floodplainA"
riverpolygonshapefile="\home\riverpolygonshapefile.shp" #a shapefile containing the "in-river" feature, i.e. the riverbed between left and right dam


#**************************************************************************
#functions block
#**************************************************************************
def vulnerabilityTotschnig(flowdepth):
    # vulnerability function (Totschnig et al. 2011)
    # Totschnig, R., Sedlacek, W., and Fuchs, S.: A quantitative vulnerability function for fluvial sediment transport
    # Nat Hazards, 58, 681-703, doi:10.1007/s11069-010-9623-5, 2011.
    # requires math
    if float(flowdepth) >=0.0:
        u = (float(flowdepth)+1.442)/1.442-1
        v = math.pow(u, 2.233)
        w = -0.443*v
        x = 1-math.exp(w)
        if x >=0.0 and x <=1.0 :
            dol = x
        elif x>1.0:
            dol = 1.0
        elif x<0.0:
            dol = 0.0
    else:
        dol = 0.0
    if dol>1:
        dol=1
    elif dol<0:
        dol=0
    return dol
def vulnerabilityPapathoma(magnitude):
    # vulnerability function (Papathoma)
    # Papathoma-Koehle, M., Zischg, A., Fuchs, S., Glade, T., and Keiler, M.:
    # Loss estimation for landslides in mountain areas - An integrated toolbox for vulnerability assessment
    # and damage documentation, Environmental Modelling & Software, 63, 156-169, doi:10.1016/j.envsoft.2014.10.003, 2015.
    # requires math
    a=-1.671
    b=3.189
    c=1.746
    if magnitude <= 0:
        dol = 0.0
    else:
        dol = 1-pow(math.e,a*pow(((magnitude+b)/b-1),c))
    if dol>1:
        dol=1.0
    elif dol<0:
        dol=0.0
    return dol
def vulnerabilityHydrotec(magnitude):
    # vulnerability function HYDROTEC (2001)
    # requires math
    if magnitude <= 0:
        dol = 0.0
    else:
        dol = (27*math.sqrt(magnitude))/100
    if dol>1:
        dol=1.0
    elif dol<0:
        dol=0.0
    return dol
def vulnerabilityJonkman(flowdepth):
    # vulnerability function JONKMAN (Jonkman 2008)
    # requires math
    vulnerabilityfile = "\home\jonkman2008.txt"
    arr=numpy.loadtxt(vulnerabilityfile,dtype='float',skiprows=1)
    depthmin=min(arr[:,0])
    depthmax=max(arr[:,0])
    dolmax=max(arr[:,1])
    x=arr[:,0]
    y=arr[:,1]
    f=interp1d(x, y)
    if flowdepth < depthmin:
        dol = 0.0
    elif flowdepth > depthmax:
        dol = dolmax
    else:
        dol = f(float(flowdepth))
    if dol>1:
        dol=1.0
    elif dol<0:
        dol=0.0
    return dol
def vulnerabilityDutta(flowdepth):
    # vulnerability function DUTTA (Dutta 2003)
    # requires math
    vulnerabilityfile = "\home\Dutta2003.txt"
    arr=numpy.loadtxt(vulnerabilityfile,dtype='float',skiprows=1)
    depthmin=min(arr[:,0])
    depthmax=max(arr[:,0])
    dolmax=max(arr[:,1])
    x=arr[:,0]
    y=arr[:,1]
    f=interp1d(x, y)
    if flowdepth < depthmin:
        dol = 0.0
    elif flowdepth > depthmax:
        dol = dolmax
    else:
        dol = f(float(flowdepth))
    if dol>1:
        dol=1.0
    elif dol<0:
        dol=0.0
    return dol
	
	
#**************************************************************************
#IBER functions
#**************************************************************************
def meshtonodesarrayIBER(meshfile):
    #reads the IBER mesh file and returns an array with NodeID, x, y, z coordinates
    mesh = open(meshfile, "r")
    #count nodes in mesh
    nodesflag = 0
    coordinatesflag=0
    countnodes = 0
    nodeslist=[]
    for line in mesh:
        checknd = line.strip().split()
        if line <> '\n' and nodesflag==0 and checknd[0] == 'MESH' and len(checknd)>=1:
            nodesflag=1
        elif nodesflag==1 and coordinatesflag == 0 and len(checknd)==1 and checknd[0]=='Coordinates':
            coordinatesflag = 1
        elif coordinatesflag == 1 and len(checknd) == 2 and checknd[0]=='end' and checknd[1]=='coordinates':
            coordinatesflag = 0
            nodesflag = 0
        elif checknd=='\n':
            coordinatesflag = 0
            nodesflag = 0
        elif nodesflag == 1 and coordinatesflag == 1 and len(checknd) == 4:
            countnodes = countnodes + 1
            nodeslist.append(checknd)
        else:
            continue
    mesh.close()
    # create an empty output array
    nodesarray = numpy.zeros((countnodes, 4), dtype=float)
    i=0
    for line in nodeslist:
        nodesarray[i,0] = line[0]
        nodesarray[i,1] = line[1]
        nodesarray[i,2] = line[2]
        nodesarray[i,3] = line[3]
        i+=1
    #check if meshnodes are in re-numbered
    if len(nodeslist)<>nodeslist[-1][0]:
        print "Attention! The mesh nodes are not numbered correctly. Please re-number the mesh"
    del (nodeslist)
    del (countnodes)
    del (line)
    del (checknd)
    return nodesarray
def createnodesshapefilefromnodesarray(nodesarray, xyztabletxtfullfilepathname, meshnodesshapefilefullpathname):
    arcpy.env.overwriteOutput = True
    #converts a nodesarray to 1. a table of NODE_ID and xyz coordinates and 2. to a point shapefile
    #numpy.savetxt(xyztabletxtfullfilepathname,nodesarray, fmt=['%1i','%1.6f','%1.6f','%1.6f'],delimiter=outputdelimiter, header="NODE_ID\tx\ty\tz")
    file=open(xyztabletxtfullfilepathname,"w")
    file.write("NODE_ID\tx\ty\tz\n")
    i=0
    for line in nodesarray:
        file.write(str(int(nodesarray[i,0]))+outputdelimiter+str(nodesarray[i,1])+outputdelimiter+str(nodesarray[i,2])+outputdelimiter+str(nodesarray[i,3])+"\n")
        i+=1
    file.close()
    outlayer = "nodesxyzLayer"
    arcpy.MakeXYEventLayer_management(table=xyztabletxtfullfilepathname, in_x_field="x", in_y_field="y", out_layer=outlayer, spatial_reference=spatialreference, in_z_field="z")
    path, file = os.path.split(meshnodesshapefilefullpathname)
    arcpy.FeatureClassToFeatureClass_conversion(outlayer, path, file)
def buildgebtonodestopologyIBER(gebshp, tempdir, meshnodesshapefilefullpathname, riverpolygonshapefile, searchradius, maximumnumberofpoints):
    #function builds a topology between the building shapefile and the nodes in nodes shapefile
    #output are an array for buildings that intersect with nodes
    #floodplainname .. name of the floodplain, used to name the output files
    #gebshp... building shapefile with obligatory fields: "Orig_FID", "Wert", "AnzPers"
    #tempdir ... directory of temporary results
    #outputdir... directory where the output will be written
    #meshfilename .msh file of mesh, full file path
    #meshnodesshapefilefullpathname .. shapefile of the nodes of the mesh with obligatory field "NODE_ID", full file path
    #riverpolygonshapefile .. shapefile of the river polygon, needed to erase those nodes located in river channel, full file path
    #method for attributing flowdepth of nodes to one building ("MAX" or "MEAN")
    #returns an array of the building shapefile, an array with the topology of intersction between nodes and buidlings, and an array of neares nodes per bulding
    #check the number of buildings
    arcpy.env.overwriteOutput = True
    countgebaeude = arcpy.GetCount_management(gebshp)
    #print "building dataset: "+gebshp+ " with " + str(countgebaeude)+ " objects"
    # **************************************************************************
    # Process: intersect buildings with nodes
    # add a field to the building shape-file and set it to 1 if the building is situated directly over >=1 nodes
    # check if column "onNode" already exists
    fieldObjList = arcpy.ListFields(gebshp)
    fieldNameList = []
    for field in fieldObjList:
        if not field.required:
            fieldNameList.append(field.name)
    if "onNode" not in fieldNameList:
        arcpy.AddField_management(gebshp, "onNode", "LONG", 2)
    arcpy.CalculateField_management(gebshp, "onNode", 0, "PYTHON")
    if "origFID" not in fieldNameList:
        arcpy.AddField_management(gebshp, "origFID", "LONG", 10)
        arcpy.CalculateField_management(gebshp, "origFID", "!FID!", "PYTHON")
    arcpy.MakeFeatureLayer_management(gebshp, 'gebaeude_lyr')
    #set onNode field to 1 if building has nodes within polygon
    #print "selecting all buildings that intersect directly with nodes ..."
    arcpy.SelectLayerByLocation_management('gebaeude_lyr', 'intersect', meshnodesshapefilefullpathname)
    arcpy.CalculateField_management('gebaeude_lyr', 'onNode', 1, "PYTHON")
    arcpy.SelectLayerByAttribute_management('gebaeude_lyr', "CLEAR_SELECTION")
    # intersect
    inFeatures = [gebshp, meshnodesshapefilefullpathname]
    intersectOutput = tempdir+"/"+"gebintersected.shp"
    #print "intersecting buildings and nodes ..."
    gebintersected = arcpy.Intersect_analysis(inFeatures, intersectOutput)
    gebintersectedarray = arcpy.da.FeatureClassToNumPyArray(gebintersected, ["origFID", "NODE_ID"])
    # **************************************************************************
    # Process: computing not-intersected buildings by near analysis
    # selects buildings that are not directly intersected
    gebnotintersected = arcpy.Select_analysis(gebshp, tempdir+"/"+"gebnotintersected.shp", '"onNode" = 0')
    #delete all nodes from nodefile that are lying within the river polygon shapefile riverpolygonshapefile
    erasenodes = tempdir+"/"+"erasenodes.shp"
    arcpy.Erase_analysis(meshnodesshapefilefullpathname, riverpolygonshapefile, erasenodes)
    # finds the 3 nearest mesh nodes around the buildings that do not intersect directly with nodes
    outtable = tempdir+"/"+"neartable.dbf"
    search_radius = str(searchradius)+' Meters'
    location = 'NO_LOCATION'
    angle = 'NO_ANGLE'
    closest = 'ALL'
    closest_count = maximumnumberofpoints
    neartable = arcpy.GenerateNearTable_analysis(gebnotintersected, erasenodes, outtable, search_radius, location, angle,closest, closest_count)
    # write original FID's and NODE_ID's to colums origFID and NODE_ID by means of joins
    arcpy.AddField_management(neartable, "origFID", "LONG", 10)
    arcpy.AddField_management(neartable, "NODE_ID", "LONG", 10)
    # make table view
    arcpy.MakeTableView_management(neartable, "neartable_view")
    # make joins, calculate fields and remove joins
    # join the gebaeude shape-file to the neartable and write the origFID into the resp. column
    arcpy.AddJoin_management("neartable_view", "IN_FID", gebnotintersected, "FID")
    calcExpression = '!gebnotintersected.origFID!'
    arcpy.CalculateField_management("neartable_view", "origFID", calcExpression, "PYTHON")
    arcpy.RemoveJoin_management("neartable_view")
    arcpy.AddJoin_management("neartable_view", "NEAR_FID", erasenodes, "FID")
    calcExpression = '"!' +'erasenodes' + '.NODE_ID!"'
    arcpy.CalculateField_management("neartable_view", "NODE_ID", calcExpression, "PYTHON")
    arcpy.RemoveJoin_management("neartable_view")
    neartablearray = arcpy.da.TableToNumPyArray(neartable, ["origFID", "NODE_ID", "NEAR_DIST"])
    #print "neartable.dbf generated ..."
    #create two lists of buildings
    #geblist = arcpy.da.FeatureClassToNumPyArray(gebshp, ["origFID", "gebCHid", "zmin50cm", "AnzPers", "Wert", "onNode"])
    # **************************************************************************
    # create and fill the array of buildings and related flow depths per timestep
    # the array has origFID of buildings-shapefile as first column and flowdepths for each timestep
    # **************************************************************************
    # search cursor for list of buildings
    gebcursor = arcpy.da.SearchCursor(gebshp, ["origFID", "gebCHid", "zmin50cm", "AnzPers", "Wert", "onNode"])
    # fill the array with origFID
    #create the gebinputarray
    # make an array from the input building table
    # array import columns: FID, Z-Min, AnzPers, Wert, area, onNode form building shapefile
    gebinputarray = numpy.zeros((int(str(countgebaeude)), 6), dtype=numpy.float32)
    j = 0
    for row in gebcursor:
        if j < int(str(countgebaeude)):
            gebinputarray[j, 0] = float(row[0])
            gebinputarray[j, 1] = float(row[1])
            gebinputarray[j, 2] = float(row[2])
            gebinputarray[j, 3] = float(row[3])
            gebinputarray[j, 4] = float(row[4])
            gebinputarray[j, 5] = float(row[5])
        j += 1
    gebcursor.reset()
    #write the files

    #delete temporary files
    arcpy.Delete_management(gebintersected)
    del(gebintersected)
    arcpy.Delete_management(intersectOutput)
    del(intersectOutput)
    arcpy.Delete_management(gebnotintersected)
    del(gebnotintersected)
    arcpy.Delete_management(erasenodes)
    del(erasenodes)
    arcpy.Delete_management(neartable)
    del(neartable)
    return gebinputarray, gebintersectedarray, neartablearray
def createflowdepthpernodeandtimesteparrayfromflowdepthfileIBER(nodesarray, iberflowdepthfile):
    #reads the flow depth results file of IBER and returns an array for each node (în first column) and timestep (colums 1 ... number of timesteps)
    flowdepthfile = open(iberflowdepthfile,"r")
    #count number of timesteps
    timestepslist=[]
    for line in flowdepthfile:
        checknd=line.strip().split()
        if len(checknd)>0 and checknd[0]=="Result" and checknd[1]=='"Depth' and checknd[2]=='(m)"' and checknd[3]=='"Hydraulic"' and checknd[-1]=='OnNodes':
            timestep=float(checknd[4])
            if timestep not in timestepslist:
                timestepslist.append(timestep)
    flowdepthfile.close()
    numberoftimesteps=len(timestepslist)
    countnodes=len(nodesarray)
    #create an array and fill it with values from flow depth results file
    timestepsarray=numpy.zeros((countnodes, numberoftimesteps+2), dtype=float) #plus one column in the beginning for NODEID and one column at the end for maximum flow depth per node
    #fill the node id's (first col)
    timestepsarray[:,0]=nodesarray[:,0]
    #create a dictionary for number of row vs NODE_ID. Is to fasten the analysis in a non-renumbered mesh
    # left is the NODE_ID and right is index
    nodesdict = {}
    i=0
    while i < len(nodesarray):
        nodesdict.update({nodesarray[i,0]: i})
        i+=1
    #fill the timestepsarray
    flowdepthfile = open(iberflowdepthfile, "r")
    depthflag=0
    valuesflag = 0
    nodeid=0
    depth=0.0
    col=1
    for line in flowdepthfile:
        depth=0.0
        nodeid = 0
        checknd=line.strip().split()
        if len(checknd)>0 and checknd[0]=="Result" and checknd[1]=='"Depth' and checknd[2]=='(m)"' and checknd[3]=='"Hydraulic"' and checknd[5]=='Scalar'and checknd[-1]=='OnNodes':
            depthflag=1
            print "timestep "+str(checknd[4])
            print col
        elif len(checknd)>0 and depthflag==1 and checknd[0]=="Values":
            valuesflag=1
        elif len(checknd)>0 and depthflag==1 and valuesflag==1 and len(checknd)==2 and checknd[0]<>"End":
            nodeid=int(checknd[0])
            depth = float(checknd[1])
            timestepsarray[int(nodesdict[nodeid]),col]=depth
        elif len(checknd)>0 and depthflag==1 and valuesflag==1 and len(checknd)==2 and checknd[0]=="End" and checknd[1]=="values":
            depthflag = 0
            valuesflag = 0
            col += 1
    flowdepthfile.close()
    #calculate last column as max of all flowdepths at timesteps
    i=0
    while i < len(timestepsarray):
        timestepsarray[i, -1] = max(timestepsarray[i, 1:-2])
        i+=1
    return timestepsarray, nodesdict, timestepslist
def calculateflowdepthperbuildingandtimestepIBER(gebinputarray, gebintersectedarray, neartablearray, timestepsarray, nodesdict, timestepslist, flowdepthstatmethod):
    #takes the array with flow depth per node and timestep, the buildings-nodes topology and computes an array for flow depths per building and timestep
    #create an array with number of rows = number of buildings and number of cols = number of timesteps +2
    gebarray=numpy.zeros((len(gebinputarray), int(len(timestepslist)) + 2), dtype=float)  # plus 1 column at the beginning for "origFID" and 1 column at the end for max flowdepth
    gebarray[:,0]=gebinputarray[:,0]
    #loop through buildings
    gebindex=0
    for geb in gebarray:
        origFID=int(gebarray[gebindex,0])
        onNode = int(gebinputarray[gebindex, 5])
        #loop through timesteps
        tsindex=1
        for timestep in timestepslist:
            flowdepthintersect=0.0
            flowdepthnotintersect=0.0
            tmpflowdepthslist=[]
            tmpflowdepthslistbuffer=[]
            tmpdistanceslist=[]
            tmpinversedistanceslist=[]
            # calculate flowdepth for buildings that have min. one node inside
            if onNode == 1:
                for itemint in gebintersectedarray:
                    if int(itemint[0]) == origFID:
                        tmpflowdepthslist.append(timestepsarray[int(nodesdict[int(itemint[1])]), tsindex])
                if flowdepthstatmethod == "MAX" and len(tmpflowdepthslist) > 0:
                    flowdepthintersect = max(tmpflowdepthslist)
                elif flowdepthstatmethod == "MEAN" and len(tmpflowdepthslist) > 0:
                    flowdepthintersect = sum(tmpflowdepthslist)/float(len(tmpflowdepthslist))
                else:
                    flowdepthintersect = -9999.0
            # calculate flowdepth for buildings on the basis of neighboring nodes
            elif onNode == 0:
                for itemnotint in neartablearray:
                    if int(itemnotint[0]) == origFID:
                        tmpflowdepthslistbuffer.append(timestepsarray[int(nodesdict[int(itemnotint[1])]), tsindex])
                        if float(itemnotint[2])==0.0:
                            tmpdistanceslist.append(0.01)
                            tmpinversedistanceslist.append(1.0/0.01)
                        else:
                            tmpdistanceslist.append(float(itemnotint[2]))
                            tmpinversedistanceslist.append(1.0/float(itemnotint[2]))
                if flowdepthstatmethod == "MAX" and len(tmpflowdepthslistbuffer) > 0:
                    flowdepthnotintersect = max(tmpflowdepthslistbuffer)
                elif flowdepthstatmethod == "MEAN" and len(tmpflowdepthslistbuffer) > 0:
                    z = 0
                    tmpflowdepth = 0.0
                    for fd in tmpflowdepthslistbuffer:
                        tmpflowdepth = tmpflowdepth + float(fd) * tmpinversedistanceslist[z]
                        z += 1
                    if sum(tmpdistanceslist) > 0:
                        flowdepthnotintersect = tmpflowdepth / sum(tmpinversedistanceslist)
                    else:
                        flowdepthnotintersect = sum(tmpflowdepthslist)/float(len(tmpflowdepthslistbuffer))
                else:
                    flowdepthnotintersect = -9999.0
            flowdepth = max(flowdepthintersect, flowdepthnotintersect)
            gebarray[gebindex,tsindex]=flowdepth
            tsindex+=1
        gebindex+=1
    #fill the last column with the maximum of flow depth of previous timesteps
    i = 0
    for row in gebarray:
        gebarray[i, -1] = max(gebarray[i, 1:-2])
        i += 1
    return gebarray
def calclossesfromgebarray(gebinputarray, gebarray):
    #creates x arrays for each vulnerability function on the beasis of the flow deth per building and timestep (gebarray) and gebinputarray
    #needs the gebinputarray and gebarray
    rows = int(gebarray.shape[0])
    cols = int(gebarray.shape[1])
    # create a damage array for the vulnerability functions
    totarr = numpy.zeros((rows, cols), dtype=numpy.float32)
    paparr = numpy.zeros((rows, cols), dtype=numpy.float32)
    hydarr = numpy.zeros((rows, cols), dtype=numpy.float32)
    jonarr = numpy.zeros((rows, cols), dtype=numpy.float32)
    dutarr = numpy.zeros((rows, cols), dtype=numpy.float32)
    totarr[:,0] = gebarray[:,0]
    paparr[:,0] = gebarray[:,0]
    hydarr[:,0] = gebarray[:,0]
    jonarr[:,0] = gebarray[:,0]
    dutarr[:,0] = gebarray[:,0]
    i = 0
    while i < rows:
        buildingvalue = float(gebinputarray[i, 4])
        j = 1
        while j < cols:
            if gebarray[i, j] > 0.0:
                totarr[i, j] = buildingvalue * vulnerabilityTotschnig(gebarray[i, j])
                paparr[i, j] = buildingvalue * vulnerabilityPapathoma(gebarray[i, j])
                hydarr[i, j] = buildingvalue * vulnerabilityHydrotec(gebarray[i, j])
                jonarr[i, j] = buildingvalue * vulnerabilityJonkman(gebarray[i, j])
                dutarr[i, j] = buildingvalue * vulnerabilityDutta(gebarray[i, j])
            j += 1
        i += 1
    # while flowdepth can decrease with the timesteps, damage cannot
    # correct decreasing values of damages
    # i.e. after peak flowdepth, the damage value would remain max
    i = 0
    j = 2
    while i < rows:
        j = 2
        while j < cols:
            if totarr[i, j] < totarr[i, j - 1]:
                totarr[i, j] = totarr[i, j - 1]
            if paparr[i, j] < paparr[i, j - 1]:
                paparr[i, j] = paparr[i, j - 1]
            if hydarr[i, j] < hydarr[i, j - 1]:
                hydarr[i, j] = hydarr[i, j - 1]
            if jonarr[i, j] < jonarr[i, j - 1]:
                jonarr[i, j] = jonarr[i, j - 1]
            if dutarr[i, j] < dutarr[i, j - 1]:
                dutarr[i, j] = dutarr[i, j - 1]
            j += 1
        i += 1
    # create a summary array
    sumarrheader = "origFID" + "\t" + "expPer" + "\t" + "damTot" + "\t" + "damPap" + "\t" + "damHyd" + "\t" + "damJon" + "\t" + "damDut" + "\t" + "damMean"
    sumarr = numpy.zeros((rows, 8), dtype=float)
    i = 0
    while i < rows:
        sumarr[i, 0] = gebarray[i, 0]
        sumarr[i, 1] = gebinputarray[i, 3]
        sumarr[i, 2] = totarr[i, -1]
        sumarr[i, 3] = paparr[i, -1]
        sumarr[i, 4] = hydarr[i, -1]
        sumarr[i, 5] = jonarr[i, -1]
        sumarr[i, 6] = dutarr[i, -1]
        sumarr[i, 7] = (totarr[i, -1] + paparr[i, -1] + hydarr[i, -1] + ngeoarr[i, -1] + jonarr[i, -1] + dutarr[i, -1]) / 5
        i += 1
    # create a total array
    totalarr = numpy.zeros((1, 8), dtype=float)
    i = 0
    anzgeb = 0
    sumpers = 0
    while i < rows:
        if sumarr[i, 7] > 0:
            anzgeb += 1
            sumpers = sumpers + sumarr[i, 1]
        i += 1
    totalarr[0, 0] = anzgeb
    totalarr[0, 1] = sumpers
    totalarr[0, 2] = sum(sumarr[:, 2])
    totalarr[0, 3] = sum(sumarr[:, 3])
    totalarr[0, 4] = sum(sumarr[:, 4])
    totalarr[0, 5] = sum(sumarr[:, 5])
    totalarr[0, 6] = sum(sumarr[:, 6])
    totalarr[0, 7] = sum(sumarr[:, 7])
    return totarr, paparr, hydarr, jonarr, dutarr, sumarr, totalarr
#**************************************************************************
#end functions block
#**************************************************************************
