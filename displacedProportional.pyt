from heapq import heappop, heappush
import arcpy, math, os, string, traceback

class Toolbox(object):
    def __init__(self):
        self.label = "Carto Visualisation"
        self.alias = "cartoVis"
        self.tools = [dorlingCartogram,lambCircular,staticProportionalSymbol,nonContigousCartogram]

class dorlingCartogram(object):
#-------------------------------------------------------------------------------
# Name:        dorling cartogram
# Purpose:
#
#   modified Python Port of original Dorling algorithm:  http://qmrg.org.uk/files/2008/11/59-area-cartograms.pdf
# and a direct python port by Zachary Forest Johnson http://indiemaps.com/code/circularCartograms/python/dorling.py
# Author:      David S. Lamb
#
# Created:     11/1/2012
# Copyright:   (c) dslamb 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
    def __init__(self):
        self.label = "Circular Cartogram (Dorling)"
        self.description = "Generates a circular cartogram based on original Dorling algorithm.  Input is polygons with a count field. http://qmrg.org.uk/files/2008/11/59-area-cartograms.pdf"
        self.friction = 0.25
        self.ratio = 0.4
        self.screen_scale = 0.001
        self.base = {}
        self.bodies = 0
        self.iters = 100
        self.structure = {}
        self.t_dist = 0.0
        self.t_radius = 0.0
        self.scaleFact = 0.0
        self.widest = 0.0

    def getParameterInfo(self):
        in_features = arcpy.Parameter(
            displayName="Input Polygons",
            name="in_features",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input")

        in_features.filter.list = ["Polygon"]

        value_field = arcpy.Parameter(
            displayName="Value Field",
            name="value_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        value_field.value = "value"
        value_field.parameterDependencies = [in_features.name]



        id_field = arcpy.Parameter(
            displayName="ID Field",
            name="id_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        id_field.value = "id"
        id_field.parameterDependencies = [in_features.name]

        iter_count = arcpy.Parameter(
            displayName = "Number of Iterations",
            name = "iter_count",
            datatype = "Long",
            parameterType = "Optional",
            direction="Input")

        iter_count.value = 100

        out_ws= arcpy.Parameter(
            displayName="Output workspace",
            name="out_workspace",
            datatype="Workspace",
            parameterType="Required",
            direction="Input")

        out_name = arcpy.Parameter(
            displayName = "Output Name",
            name ="out_name",
            datatype = "String",
            parameterType="Required",
            direction="Input"
        )


        parameters =[in_features, value_field, id_field, iter_count, out_ws, out_name]
        return parameters

    def updateParameters(self, params):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parmater
        has been changed."""

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        if parameters[0].value:
            fc= parameters[0].value
            desc = arcpy.Describe(fc)
            if hasattr(desc, "spatialReference"):
                if desc.spatialReference.PCSCode == 0:
                    parameters[0].setErrorMessage("Should be a projected coordinate system (not latitude and longitude).")

    def execute(self, parameters, messages):
            lyr = parameters[0].value
            arcpy.env.workspace = parameters[4].value
            arcpy.env.overwriteOutput = True
            flds = arcpy.ListFields(lyr)
            fld = arcpy.Field()
            fldID = arcpy.Field()
            done = 0
            for f in flds:
                arcpy.AddMessage(f.name)
                if (f.name == parameters[1].valueAsText):
                    arcpy.AddMessage(parameters[1].value)
                    fld = f
                if(f.name == parameters[2].valueAsText):
                    fldID = f

            self.bodies = int(arcpy.GetCount_management(lyr).getOutput(0))
            arcpy.AddMessage("Bodies %s" %(self.bodies))
            self.iters = parameters[3].value
            desc = arcpy.Describe(lyr)
            clause = arcpy.AddFieldDelimiters(lyr, fld.name) + " > 0"
            cursor = arcpy.da.SearchCursor(lyr, ["SHAPE@", "SHAPE@XY", parameters[1].valueAsText, parameters[2].valueAsText],where_clause=clause)
            if (cursor):
                arcpy.AddMessage(fld)
                outName =  parameters[5].valueAsText
                fc = arcpy.CreateFeatureclass_management(arcpy.env.workspace, outName, "POLYGON",spatial_reference = desc.spatialReference)
                arcpy.AddField_management(fc, fld.name, fld.type, fld.precision, fld.scale, fld.length)
                arcpy.AddField_management(fc, fldID.name, fldID.type, fldID.precision, fldID.scale, fldID.length)
                keyID=0
                for row in cursor:
                    self.base[keyID] = [row[0],row[1][0],row[1][1], row[2], row[3]]
                    keyID += 1
                del cursor
                self.buildRelationsList(messages)
                self.scaleFact = self.t_dist / self.t_radius
                arcpy.AddMessage("tdist %s tradius %s" % (self.t_dist, self.t_radius))

                for k in self.structure.iterkeys():
                    self.structure[k].radius = self.scaleFact * math.sqrt(self.structure[k].value/math.pi)
                    if self.structure[k].radius > self.widest:
                        self.widest = self.structure[k].radius

                #arcpy.AddMessage("Scale Factor: %s, Widest: %s" %(self.scaleFact, self.widest))
                for i in range(self.iters):
                    displacement = 0.0
                    for k in self.structure.iterkeys():
                        number = 0
                        struct = self.structure[k]
                        distance = self.widest + struct.radius
                        xrepel = 0.0
                        yrepel = 0.0
                        xattract = 0.0
                        yattract = 0.0
                        xtotal = 0.0
                        ytotal = 0.0
                        closest = self.widest

                        #k_buff = struct.getArcPointGeometry().buffer(struct.radius)
                        #arcpy.AddMessage("Current k: %s radius %s" %(k, struct.radius))
                        i = 0
                        for j in self.structure.iterkeys():
                            if k != j:
                                j_struct =self.structure[j]
                                #arcpy.AddMessage("k: %s j: %s" %(k,j))
                                if (struct.new_x-distance < j_struct.new_x and	struct.new_x+distance >= j_struct.new_x):
                                    if (struct.new_y-distance < j_struct.new_y and struct.new_y+distance >= j_struct.new_y):
                                        i+=1
                                        kpnt = struct.getArcPointGeometry()
                                        jpnt = self.structure[j].getArcPointGeometry()
                                        dist = kpnt.distanceTo(jpnt)
                                        if dist < closest:
                                            closest = dist
                                        overlap = struct.radius + j_struct.radius - dist
                                        if overlap > 0.0:
                                            if dist > 1.0:
                                                xrepel = xrepel - overlap*(j_struct.new_x-struct.new_x)/dist
                                                yrepel = yrepel - overlap*(j_struct.new_y - struct.new_y)/dist

                        for j in struct.boundaries.iterkeys():
                            j_struct =self.structure[j]
                            dist = struct.getArcPointGeometry().distanceTo(j_struct.getArcPointGeometry())
                            overlapAttract = dist - struct.radius - j_struct.radius
                            if overlapAttract > 0.0:
                                overlapAttract = overlapAttract * struct.boundaries[j]/struct.perimeter
                                xattract = xattract + overlapAttract * (j_struct.new_x - struct.new_x)/dist
                                yattract = yattract + overlapAttract * (j_struct.new_y - struct.new_y)/dist

                        atrdst = math.sqrt(xattract * xattract + yattract * yattract)
                        repdst = math.sqrt(xrepel * xrepel + yrepel * yrepel)
                        if repdst > closest:
                            xrepel = closest * xrepel / (repdst + 1.0)
                            yrepel = closest * yrepel / (repdst + 1.0)
                            repdst = closest
                        if repdst > 0.0:
                            xtotal = (1.0-self.ratio) * xrepel + self.ratio * (repdst * xattract / (atrdst + 1.0))
                            ytotal = (1.0-self.ratio) * yrepel + self.ratio * (repdst * yattract / (atrdst + 1.0))
                        else:
                            if atrdst > closest:
                                xattract = closest * xattract/(atrdst+1.0)
                                yattract = closest * yattract/(atrdst+1.0)
                            xtotal = xattract
                            ytotal = yattract

                        struct.disp_x = self.friction * (struct.disp_x + xtotal)
                        struct.disp_y = self.friction * (struct.disp_y + ytotal)

                    for k in self.structure.iterkeys():
                        struct = self.structure[k]
                        struct.new_x += struct.disp_x + 0.5
                        struct.new_y += struct.disp_y + 0.5

                    #done += 1
                    #displacement = displacement / bodies


                i_cur = arcpy.da.InsertCursor(fc, ["SHAPE@", fld.name, fldID.name])
                for k in self.structure.iterkeys():
                    buff = self.structure[k].getArcPointGeometry().buffer(self.structure[k].radius)
                    i_cur.insertRow([buff, self.structure[k].value, self.structure[k].ident])
                del i_cur



    def buildRelationsList(self, messages):
        try:
            arcpy.AddMessage("Build relations list")
            for k in self.base.iterkeys():
                k_geom = self.base[k][0]
                if not self.structure.has_key(k):
                    struct = relationObject(self.base[k][4],int(self.base[k][3]))
                    struct.orig_x = int(self.base[k][1])
                    struct.orig_y = int(self.base[k][2])
                    struct.new_x = struct.orig_x
                    struct.new_y = struct.orig_y
                    self.structure[k] = struct
                    #self.structure[k].perimeter = k_geom.length
                for j in self.base.iterkeys():
                    if k != j:
                        if k_geom.touches(self.base[j][0]):
                        ##arcpy.AddMessage("%s Touches %s" %(k,j))
                            inter_geom = k_geom.intersect(self.base[j][0], 2)
                            self.structure[k].addBoundary(j, inter_geom.length)
                            self.structure[k].perimeter += inter_geom.length
                        ##arcpy.AddMessage( "%s touches %s with length %s" %(k,j,inter_geom.length))
            i = 1
            for k in self.structure.iterkeys():
                for j in self.structure[k].boundaries.iterkeys():
                    if j>0:
                        if j < k:
                            xd = self.structure[k].orig_x - self.structure[j].orig_x
                            yd = self.structure[k].orig_y - self.structure[j].orig_y
                            self.t_dist += math.sqrt(xd*xd+yd*yd)
                            self.t_radius += math.sqrt(self.structure[k].value/math.pi) + math.sqrt(self.structure[j].value/math.pi)
                            i += 1
        except:
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError("While building a relationship list the following errors occured:")
            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)





class lambCircular(object):
    #-------------------------------------------------------------------------------
# Name:        displaced proportional symbol map for circles
# Purpose:
#
#   sources: http://arxiv.org/pdf/0911.0626.pdf
#http://www.isprs.org/proceedings/XXXVIII/part4/files/Inoue.pdf
#http://www.graphviz.org/Documentation/GH10.pdf
#http://research.microsoft.com/en-us/um/people/weiweicu/images/energybundling.pdf
#http://www.cse.ust.hk/~zhouhong/articles/infovis08_weiwei.pdf
#http://www.tandfonline.com/doi/pdf/10.1080/13658810802186882
# Author:      David S Lamb
#
# Created:     11/1/2012
# Copyright:   (c) dslamb 2012
# Licence:     <your licence>
#-----------------------
    def __init__(self):
        self.label = "Displaced Proportional Symbol, Circles (PRISM)"
        self.description = "Generates proportional symbols, and removes overlaps using PRoxImity Stress Model (PRISM) algorithm for layout. http://arxiv.org/pdf/0911.0626.pdf"

    def updateParameters(self, params):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parmater
        has been changed."""

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        if parameters[0].value:
            fc= parameters[0].value
            desc = arcpy.Describe(fc)
            if hasattr(desc, "spatialReference"):
                if desc.spatialReference.PCSCode == 0:
                    parameters[0].setErrorMessage("Should be a projected coordinate system (not latitude and longitude).")

    def getParameterInfo(self):
        in_features = arcpy.Parameter(
            displayName="Input Features",
            name="in_features",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input")

        in_features.filter.list = ["Polygon", "Point"]

        value_field = arcpy.Parameter(
            displayName="Value Field (Zero values indicate missing data)",
            name="value_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        value_field.value = "value"
        value_field.parameterDependencies = [in_features.name]

        id_field = arcpy.Parameter(
            displayName="ID Field",
            name="id_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        id_field.value = "id"
        id_field.parameterDependencies = [in_features.name]

        symb_size= arcpy.Parameter(
            displayName="Minimum circle radius in coordinate system units (feet, meters, etc..)",
            name="symb_size",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        symb_size.value = "100"

        out_ws= arcpy.Parameter(
            displayName="Output workspace",
            name="out_workspace",
            datatype="Workspace",
            parameterType="Required",
            direction="Input")

        out_fn= arcpy.Parameter(
            displayName="Output Filename",
            name="out_fn",
            datatype="String",
            parameterType="Required",
            direction="Input")

        parameters =[in_features, value_field, id_field, symb_size,  out_ws, out_fn]
        return parameters

    def execute(self, parameters, messages):
        try:
            inFeatures  = parameters[0].valueAsText
            #Name of Value Field
            valnm = parameters[1].valueAsText
            #Name of ID Field
            idnm = parameters[2].valueAsText

            min_symbol_size = parameters[3].value

            #Output workspace
            outFeatures   = parameters[4].valueAsText
            #Output Name
            outName   = parameters[5].valueAsText

            arcpy.env.overwriteOutput = True
            desc = arcpy.Describe(inFeatures)

            #hold the actual fields
            valFld = arcpy.Field()
            idFld = arcpy.Field()
            for fld in arcpy.ListFields(inFeatures):
                if fld.name == valnm:
                    valFld = fld
                if fld.name == idnm:
                    idFld = fld


            if valFld and idFld:
                arcpy.AddMessage("Fields Exist")

                #Create the empty output field and add the field information
                outFC = arcpy.CreateFeatureclass_management(outFeatures,outName, "POLYGON", spatial_reference=desc.spatialReference)
                arcpy.AddField_management(outFC, valFld.name, valFld.type, valFld.precision, valFld.scale, valFld.length)
                arcpy.AddField_management(outFC, idFld.name, idFld.type, idFld.precision, idFld.scale, idFld.length)

                #Create a featureset for speed
                inFS = arcpy.FeatureSet(inFeatures)
                #Can't handle values of zero
                clause = arcpy.AddFieldDelimiters(inFeatures, valFld.name) + " > 0"
                max_value = 0.0
                min_value = 1000000000^10
                max_sqrt_value = 0.0
                min_sqrt_value = 0.0
                t_dist = 0
                t_radius = 0
                with arcpy.da.SearchCursor(inFS, ["SHAPE@XY", idFld.name, valFld.name],where_clause=clause) as sc:
                    nodes = []
                    i = 0

                    for row in sc:
                        nd = Node(i, row[2], row[1])
                        max_value = max(max_value,row[2])
                        if row[2] >0:
                            min_value = min(min_value,row[2])
                        nd.x = row[0][0]
                        nd.y = row[0][1]
                        nodes.append(nd)
                        i+=1

                arcpy.AddMessage("Minimum Value: %s"%(min_value))
                if max_value > 0:
                    max_sqrt_value = math.sqrt(max_value)
                if min_value > 0:
                    min_sqrt_value = math.sqrt(min_value)
                else:
                    arcpy.AddError("No values over 0.0")
                arcpy.AddMessage("Minimum Square Root Value: %s"%(min_sqrt_value))
                #scale = max((desc.extent.XMax-desc.extent.XMin),(desc.extent.YMax-desc.extent.YMin)) / t_radius
                #for n in nodes:
                    #n.radius *= scale

                for nd in nodes:
                    nd.sqrt_value = math.sqrt(nd.value)
                    nd.radius = min_symbol_size * float(nd.sqrt_value)/min_sqrt_value
                    nd.area = math.pi * (nd.radius*nd.radius)

                nodes.sort(key=lambda x: x.area, reverse=True)
                arcpy.AddMessage(str(len(nodes)))
                testPSM = True
                while testPSM:
                    G,tri = utils.buildDelaunayGraph(nodes)
                    d = G.getNodes()
                    psm = 0.0
                    for n in d.iterkeys():
                        node = d[n]

                        for k in node.getNeighbors().iterkeys():
                            other = G.getNode(k)
                            dist = max(node.getNeighbors()[k], 0.000001)
                            #fac = float(node.radius + other.radius)/(dist)
                            try:
                                fac = float(node.radius+other.radius)/dist#(node.radius + other.radius)/dist
                                hfac = float(node.radius+other.radius)/dist
                            except ZeroDivisionError:
                                fac = float(node.radius+other.radius)/.0000001#(node.radius + other.radius)/dist
                                hfac = float(node.radius+other.radius)/.0000001


                            fac = min(fac,hfac)
                            tij = max(fac, 1)
                            smax = 1.5
                            sij = min(tij, smax)
                            dij = sij*dist
                            if dist >0:
                                newX = node.x + (other.x-node.x) / dist * dij
                                newY = node.y + (other.y-node.y) / dist * dij
                                other.disp_x +=newX - other.x
                                other.disp_y += newY - other.y

                            psm += (1/(dist*dist)) * ((dist-dij)*(dist-dij))
                    arcpy.AddMessage("Current psm %s"%(psm))
                    nodes = []
                    for n in d.iterkeys():
                        node = d[n]
                        node.x += node.disp_x
                        node.y += node.disp_y
                        node.disp_x = 0.0
                        node.disp_y = 0.0
                        node.clearNeighbors()
                        nodes.append(node)
                    if psm >= 0.000001:
                        testPSM = True
                    else:
                        testPSM = False
                #check for overlaps amongst everything
                testPSM = True
                while testPSM:
                    G, tri = utils.buildDelaunayGraph(nodes)
                    d = G.getNodes()
                    psm = 0.0
                    for n in d.iterkeys():
                        node = d[n]
                        for k in d.iterkeys():
                            if n!=k:
                                other = G.getNode(k)
                                dist = max(utils.lengthBetweenNodes(node,other), 0.000001)
                                #fac = float(node.radius + other.radius)/(dist)#(node.radius + other.radius)/dist
                                try:
                                    fac = float(node.radius+other.radius)/dist#(node.radius + other.radius)/dist
                                    hfac = float(node.radius+other.radius)/dist
                                except ZeroDivisionError:
                                    fac = float(node.radius+other.radius)/.0000001#(node.radius + other.radius)/dist
                                    hfac = float(node.radius+other.radius)/.0000001
                                fac = min(fac,hfac)
                                tij = max(fac, 1)
                                smax = 1.5
                                sij = min(tij, smax)
                                dij = sij*dist
                                if dist >0:
                                    newX = node.x + (other.x-node.x) / dist * dij
                                    newY = node.y + (other.y-node.y) / dist * dij
                                    other.disp_x +=newX - other.x
                                    other.disp_y += newY - other.y

                                psm += (1/(dist*dist)) * ((dist-dij)*(dist-dij))
                    arcpy.AddMessage("Second loop psm %s"%(psm))
                    nodes = []
                    for n in d.iterkeys():
                        node = d[n]
                        node.x += node.disp_x
                        node.y += node.disp_y
                        node.disp_x = 0.0
                        node.disp_y = 0.0
                        node.clearNeighbors()
                        nodes.append(node)
                    if psm >= 0.000001:
                        testPSM = True
                    else:
                        testPSM = False
                ic = arcpy.da.InsertCursor(outFC, ["SHAPE@",idFld.name, valFld.name])

                for n in nodes:
                    pnt = arcpy.Point(n.x, n.y)
                    buf = arcpy.PointGeometry(pnt).buffer(n.radius)
                    #arcpy.AddMessage(buf)
                    ic.insertRow([buf, n.altid, n.value])
                del ic
        except:
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError("While creating a circular cartogram the following errors occured:")
            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)

class lambSquare(object):
    #-------------------------------------------------------------------------------
# Name:        displaced proportional symbol map for squares
# Purpose:
#
#   sources: http://arxiv.org/pdf/0911.0626.pdf
#http://www.isprs.org/proceedings/XXXVIII/part4/files/Inoue.pdf
#http://www.graphviz.org/Documentation/GH10.pdf
#http://research.microsoft.com/en-us/um/people/weiweicu/images/energybundling.pdf
#http://www.cse.ust.hk/~zhouhong/articles/infovis08_weiwei.pdf
#http://www.tandfonline.com/doi/pdf/10.1080/13658810802186882
# Author:      David S Lamb
#
# Created:     11/1/2012
# Copyright:   (c) dslamb 2012
# Licence:     <your licence>
#-----------------------
    def __init__(self):
        self.label = "Displaced Proportional Symbol, Squares (PRISM)"
        self.description = "Generates proportional symbols, and removes overlaps using PRoxImity Stress Model (PRISM) algorithm for layout. http://arxiv.org/pdf/0911.0626.pdf"

    def updateParameters(self, params):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parmater
        has been changed."""

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        if parameters[0].value:
            fc= parameters[0].value
            desc = arcpy.Describe(fc)
            if hasattr(desc, "spatialReference"):
                if desc.spatialReference.PCSCode == 0:
                    parameters[0].setErrorMessage("Should be a projected coordinate system (not latitude and longitude).")

    def getParameterInfo(self):
        in_features = arcpy.Parameter(
            displayName="Input Features",
            name="in_features",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input")

        in_features.filter.list = ["Polygon", "Point"]

        value_field = arcpy.Parameter(
            displayName="Value Field (Zero values indicate missing data)",
            name="value_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        value_field.value = "value"
        value_field.parameterDependencies = [in_features.name]

        id_field = arcpy.Parameter(
            displayName="ID Field",
            name="id_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        id_field.value = "id"
        id_field.parameterDependencies = [in_features.name]

        symb_size= arcpy.Parameter(
            displayName="Minimum circle radius in coordinate system units (feet, meters, etc..)",
            name="symb_size",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        symb_size.value = "100"

        out_ws= arcpy.Parameter(
            displayName="Output workspace",
            name="out_workspace",
            datatype="Workspace",
            parameterType="Required",
            direction="Input")

        out_fn= arcpy.Parameter(
            displayName="Output Filename",
            name="out_fn",
            datatype="String",
            parameterType="Required",
            direction="Input")

        parameters =[in_features, value_field, id_field, symb_size,  out_ws, out_fn]
        return parameters

    def execute(self, parameters, messages):
        try:
            inFeatures  = parameters[0].valueAsText
            #Name of Value Field
            valnm = parameters[1].valueAsText
            #Name of ID Field
            idnm = parameters[2].valueAsText

            min_symbol_size = parameters[3].value

            #Output workspace
            outFeatures   = parameters[4].valueAsText
            #Output Name
            outName   = parameters[5].valueAsText

            arcpy.env.overwriteOutput = True
            desc = arcpy.Describe(inFeatures)

            #hold the actual fields
            valFld = arcpy.Field()
            idFld = arcpy.Field()
            for fld in arcpy.ListFields(inFeatures):
                if fld.name == valnm:
                    valFld = fld
                if fld.name == idnm:
                    idFld = fld


            if valFld and idFld:
                arcpy.AddMessage("Fields Exist")

                #Create the empty output field and add the field information
                outFC = arcpy.CreateFeatureclass_management(outFeatures,outName, "POLYGON", spatial_reference=desc.spatialReference)
                arcpy.AddField_management(outFC, valFld.name, valFld.type, valFld.precision, valFld.scale, valFld.length)
                arcpy.AddField_management(outFC, idFld.name, idFld.type, idFld.precision, idFld.scale, idFld.length)

                #Create a featureset for speed
                inFS = arcpy.FeatureSet(inFeatures)
                #Can't handle values of zero
                clause = arcpy.AddFieldDelimiters(inFeatures, valFld.name) + " > 0"
                max_value = 0.0
                min_value = 1000000000^10
                max_sqrt_value = 0.0
                min_sqrt_value = 0.0
                t_dist = 0
                t_radius = 0
                with arcpy.da.SearchCursor(inFS, ["SHAPE@XY", idFld.name, valFld.name],where_clause=clause) as sc:
                    nodes = []
                    i = 0

                    for row in sc:
                        nd = Node(i, row[2], row[1])
                        max_value = max(max_value,row[2])
                        if row[2] >0:
                            min_value = min(min_value,row[2])
                        nd.x = row[0][0]
                        nd.y = row[0][1]
                        nodes.append(nd)
                        i+=1

                arcpy.AddMessage("Minimum Value: %s"%(min_value))
                if max_value > 0:
                    max_sqrt_value = math.sqrt(max_value)
                if min_value > 0:
                    min_sqrt_value = math.sqrt(min_value)
                else:
                    arcpy.AddError("No values over 0.0")
                arcpy.AddMessage("Minimum Square Root Value: %s"%(min_sqrt_value))
                #scale = max((desc.extent.XMax-desc.extent.XMin),(desc.extent.YMax-desc.extent.YMin)) / t_radius
                #for n in nodes:
                    #n.radius *= scale

                for nd in nodes:
                    nd.sqrt_value = math.sqrt(nd.value)
                    nd.radius = min_symbol_size * float(nd.sqrt_value)/min_sqrt_value


                arcpy.AddMessage(str(len(nodes)))
                testPSM = True
                while testPSM:
                    G,tri = utils.buildDelaunayGraph(nodes)
                    d = G.getNodes()
                    psm = 0.0
                    for n in d.iterkeys():
                        node = d[n]

                        for k in node.getNeighbors().iterkeys():
                            other = G.getNode(k)
                            dist = max(node.getNeighbors()[k], 0.000001)
                            #fac = float(node.radius + other.radius)/(dist)
                            try:
                                fac = float(node.radius+other.radius)/dist#(node.radius + other.radius)/dist
                                hfac = float(node.radius+other.radius)/dist
                            except ZeroDivisionError:
                                fac = float(node.radius+other.radius)/.0000001#(node.radius + other.radius)/dist
                                hfac = float(node.radius+other.radius)/.0000001


                            fac = min(fac,hfac)
                            tij = max(fac, 1)
                            smax = 1.5
                            sij = min(tij, smax)
                            dij = sij*dist
                            if dist >0:
                                newX = node.x + (other.x-node.x) / dist * dij
                                newY = node.y + (other.y-node.y) / dist * dij
                                other.disp_x +=newX - other.x
                                other.disp_y += newY - other.y

                            psm += (1/(dist*dist)) * ((dist-dij)*(dist-dij))
                    arcpy.AddMessage("Current psm %s"%(psm))
                    nodes = []
                    for n in d.iterkeys():
                        node = d[n]
                        node.x += node.disp_x
                        node.y += node.disp_y
                        node.disp_x = 0.0
                        node.disp_y = 0.0
                        node.clearNeighbors()
                        nodes.append(node)
                    if psm >= 0.000001:
                        testPSM = True
                    else:
                        testPSM = False
                #check for overlaps amongst everything
                testPSM = True
                while testPSM:
                    G, tri = utils.buildDelaunayGraph(nodes)
                    d = G.getNodes()
                    psm = 0.0
                    for n in d.iterkeys():
                        node = d[n]
                        for k in d.iterkeys():
                            if n!=k:
                                other = G.getNode(k)
                                dist = max(utils.lengthBetweenNodes(node,other), 0.000001)
                                #fac = float(node.radius + other.radius)/(dist)#(node.radius + other.radius)/dist
                                try:
                                    fac = float(node.radius+other.radius)/dist#(node.radius + other.radius)/dist
                                    hfac = float(node.radius+other.radius)/dist
                                except ZeroDivisionError:
                                    fac = float(node.radius+other.radius)/.0000001#(node.radius + other.radius)/dist
                                    hfac = float(node.radius+other.radius)/.0000001
                                fac = min(fac,hfac)
                                tij = max(fac, 1)
                                smax = 1.5
                                sij = min(tij, smax)
                                dij = sij*dist
                                if dist >0:
                                    newX = node.x + (other.x-node.x) / dist * dij
                                    newY = node.y + (other.y-node.y) / dist * dij
                                    other.disp_x +=newX - other.x
                                    other.disp_y += newY - other.y

                                psm += (1/(dist*dist)) * ((dist-dij)*(dist-dij))
                    arcpy.AddMessage("Second loop psm %s"%(psm))
                    nodes = []
                    for n in d.iterkeys():
                        node = d[n]
                        node.x += node.disp_x
                        node.y += node.disp_y
                        node.disp_x = 0.0
                        node.disp_y = 0.0
                        node.clearNeighbors()
                        nodes.append(node)
                    if psm >= 0.000001:
                        testPSM = True
                    else:
                        testPSM = False
                ic = arcpy.da.InsertCursor(outFC, ["SHAPE@",idFld.name, valFld.name])

                for n in nodes:
                    pnt = arcpy.Point(n.x, n.y)
                    buf = arcpy.PointGeometry(pnt).buffer(n.radius)
                    #arcpy.AddMessage(buf)
                    ic.insertRow([buf, n.altid, n.value])
                del ic
        except:
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError("While creating a circular cartogram the following errors occured:")
            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)



class staticProportionalSymbol(object):
    #-------------------------------------------------------------------------------
# Name:        Static Proportional Symbol
# Purpose:
# Author:      David S Lamb
#
# Created:     11/1/2012
# Copyright:   (c) dslamb 2012
# Licence:     <your licence>
#-----------------------
    def __init__(self):
        self.label = "Static Proportional Symbol (Circle)"
        self.description = "Generates proportional symbols"

    def updateParameters(self, params):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parmater
        has been changed."""

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        if parameters[0].value:
            fc= parameters[0].value
            desc = arcpy.Describe(fc)
            if hasattr(desc, "spatialReference"):
                if desc.spatialReference.PCSCode == 0:
                    parameters[0].setErrorMessage("Should be a projected coordinate system (not latitude and longitude).")

    def getParameterInfo(self):
        in_features = arcpy.Parameter(
            displayName="Input Features",
            name="in_features",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input")

        in_features.filter.list = ["Polygon", "Point"]

        value_field = arcpy.Parameter(
            displayName="Value Field (Zero values indicate missing data)",
            name="value_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        value_field.value = "value"
        value_field.parameterDependencies = [in_features.name]

        id_field = arcpy.Parameter(
            displayName="ID Field",
            name="id_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        id_field.value = "id"
        id_field.parameterDependencies = [in_features.name]

        symb_size= arcpy.Parameter(
            displayName="Minimum circle radius in coordinate system units (feet, meters, etc..)",
            name="symb_size",
            datatype="Long",
            parameterType="Required",
            direction="Input")
        symb_size.value = "100"

        out_ws= arcpy.Parameter(
            displayName="Output workspace",
            name="out_workspace",
            datatype="Workspace",
            parameterType="Required",
            direction="Input")

        out_fn= arcpy.Parameter(
            displayName="Output Filename",
            name="out_fn",
            datatype="String",
            parameterType="Required",
            direction="Input")

        parameters =[in_features, value_field, id_field, symb_size,  out_ws, out_fn]
        return parameters

    def execute(self, parameters, messages):
        try:
            inFeatures  = parameters[0].valueAsText
            #Name of Value Field
            valnm = parameters[1].valueAsText
            #Name of ID Field
            idnm = parameters[2].valueAsText

            min_symbol_size = parameters[3].value

            #Output workspace
            outFeatures   = parameters[4].valueAsText
            #Output Name
            outName   = parameters[5].valueAsText

            arcpy.env.overwriteOutput = True
            desc = arcpy.Describe(inFeatures)

            #hold the actual fields
            valFld = arcpy.Field()
            idFld = arcpy.Field()
            for fld in arcpy.ListFields(inFeatures):
                if fld.name == valnm:
                    valFld = fld
                if fld.name == idnm:
                    idFld = fld


            if valFld and idFld:
                arcpy.AddMessage("Fields Exist")

                #Create the empty output field and add the field information
                outFC = arcpy.CreateFeatureclass_management(outFeatures,outName, "POLYGON", spatial_reference=desc.spatialReference)
                arcpy.AddField_management(outFC, valFld.name, valFld.type, valFld.precision, valFld.scale, valFld.length)
                arcpy.AddField_management(outFC, idFld.name, idFld.type, idFld.precision, idFld.scale, idFld.length)

                #Create a featureset for speed
                inFS = arcpy.FeatureSet(inFeatures)
                #Can't handle values of zero
                clause = arcpy.AddFieldDelimiters(inFeatures, valFld.name) + " > 0"
                max_value = 0.0
                min_value = 1000000000^10
                max_sqrt_value = 0.0
                min_sqrt_value = 0.0
                t_dist = 0
                t_radius = 0
                with arcpy.da.SearchCursor(inFS, ["SHAPE@XY", idFld.name, valFld.name],where_clause=clause) as sc:
                    nodes = []
                    i = 0

                    for row in sc:
                        nd = Node(i, row[2], row[1])
                        max_value = max(max_value,row[2])
                        if row[2] >0:
                            min_value = min(min_value,row[2])
                        nd.x = row[0][0]
                        nd.y = row[0][1]
                        nodes.append(nd)
                        i+=1

                arcpy.AddMessage("Minimum Value: %s"%(min_value))
                if max_value > 0:
                    max_sqrt_value = math.sqrt(max_value)
                if min_value > 0:
                    min_sqrt_value = math.sqrt(min_value)
                else:
                    arcpy.AddError("No values over 0.0")
                arcpy.AddMessage("Minimum Square Root Value: %s"%(min_sqrt_value))
                #scale = max((desc.extent.XMax-desc.extent.XMin),(desc.extent.YMax-desc.extent.YMin)) / t_radius
                #for n in nodes:
                    #n.radius *= scale

                for nd in nodes:
                    nd.sqrt_value = math.sqrt(nd.value)
                    nd.radius = min_symbol_size * float(nd.sqrt_value)/min_sqrt_value
                    nd.area = math.pi * (nd.radius*nd.radius)


                ic = arcpy.da.InsertCursor(outFC, ["SHAPE@",idFld.name, valFld.name])
                #sort nodes by area of circle, descending order
                nodes.sort(key=lambda x: x.area, reverse=True)
                for n in nodes:
                    pnt = arcpy.Point(n.x, n.y)
                    buf = arcpy.PointGeometry(pnt).buffer(n.radius)
                    #arcpy.AddMessage(buf)
                    ic.insertRow([buf, n.altid, n.value])
                del ic
        except:
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError("While creating a circular cartogram the following errors occured:")
            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)

class nonContigousCartogram(object):
    #-------------------------------------------------------------------------------
# Name:        nonContiguous cartogram
# Purpose:  Creates a nonContiguous cartogram
#
# Author:      David S Lamb
#
# Created:     11/1/2012
# Copyright:   (c) dslamb 2012
# Licence:     <your licence>
#-----------------------
    def __init__(self):
        self.label = "Noncontiguous Cartogram"
        self.description = "Generates a noncontiguous cartogram (scales polygons according to the value), then uses PRoxImity Stress Model (PRISM) algorithm for layout. http://arxiv.org/pdf/0911.0626.pdf"
        self.count = 0
    def updateParameters(self, params):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parmater
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        if parameters[0].value:
            fc= parameters[0].value
            desc = arcpy.Describe(fc)
            if hasattr(desc, "spatialReference"):
                if desc.spatialReference.PCSCode == 0:
                    parameters[0].setErrorMessage("Should be a projected coordinate system (not latitude and longitude).")
        return
    def getParameterInfo(self):
        in_features = arcpy.Parameter(
            displayName="Input Features",
            name="in_features",
            datatype="Feature Layer",
            parameterType="Required",
            direction="Input")

        in_features.filter.list = ["Polygon"]

        value_field = arcpy.Parameter(
            displayName="Value Field",
            name="value_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        value_field.parameterDependencies = [in_features.name]

        id_field = arcpy.Parameter(
            displayName="ID Field",
            name="id_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")

        id_field.parameterDependencies = [in_features.name]


        out_ws= arcpy.Parameter(
            displayName="Output workspace",
            name="out_workspace",
            datatype="Workspace",
            parameterType="Required",
            direction="Input")

        out_fn= arcpy.Parameter(
            displayName="Output Filename",
            name="out_fn",
            datatype="String",
            parameterType="Required",
            direction="Input")

        outline_bool = arcpy.Parameter(
            displayName="Use feature extent instead of shape",
            name="outline_bool",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")

        parameters =[in_features, value_field, id_field,  out_ws, out_fn,outline_bool]
        return parameters

    def execute(self, parameters, messages):
        try:
            inFeatures  = parameters[0].valueAsText
            valnm = parameters[1].valueAsText
            idnm = parameters[2].valueAsText
            outFeatures   = parameters[3].valueAsText
            outName   = parameters[4].valueAsText
            outline = parameters[5].value

            arcpy.env.overwriteOutput = True
            desc = arcpy.Describe(inFeatures)
            valFld = None
            idFld = None
            for fld in arcpy.ListFields(inFeatures):
                if fld.name == valnm:
                    valFld = fld
                if fld.name == idnm:
                    idFld = fld
            if valFld and idFld:
                arcpy.AddMessage("Fields Exist")
                arcpy.AddMessage(idFld.name)
                outFC = arcpy.CreateFeatureclass_management(outFeatures,outName, "POLYGON", spatial_reference=desc.spatialReference)
                arcpy.AddField_management(outFC, valFld.name, valFld.type, valFld.precision, valFld.scale, valFld.length)
                arcpy.AddField_management(outFC, idFld.name, idFld.type, idFld.precision, idFld.scale, idFld.length)

                inFS = arcpy.FeatureSet(inFeatures)
                clause = arcpy.AddFieldDelimiters(inFeatures, valFld.name) + " > 0"
                sc = arcpy.da.SearchCursor(inFS, ["SHAPE@XY","SHAPE@", idFld.name, valFld.name],where_clause=clause)
                nodes = []
                i = 0
                maxValue = 0.0
                maxArea = 0.0
                minArea = 0.0
                minValue = math.pow(10,10)
                for row in sc:
                    nd = Node(i, row[3], row[2])
                    if nd.value > maxValue:
                        maxValue = nd.value
                    if nd.value < minValue:
                        minValue = nd.value

                    nd.x = row[0][0]
                    nd.y = row[0][1]
                    nd.o_x = nd.x
                    nd.o_y = nd.y
                    nd.geometry = row[1]
                    if nd.geometry.area > maxArea:
                        maxArea = nd.geometry.area

                    nodes.append(nd)
                    i+=1

                minArea = (float(minValue) / maxValue) * maxArea
                arcpy.AddMessage("Resizing Geometry...")
                for n in nodes:
                    n.radius = (float(n.value) / maxValue) * maxArea
                    n.scale_fact = math.sqrt(float(n.radius)/n.geometry.area)
                    n.geometry = utils.resizeGeometry(n.geometry, n.scale_fact, minArea)
                    if n.geometry:
                        n.x = n.o_x = n.geometry.centroid.X
                        n.y = n.o_y = n.geometry.centroid.Y
                        n.width = n.geometry.extent.XMax - n.x
                        n.height = n.geometry.extent.YMax - n.y

                arcpy.AddMessage("Improving Layout...")
                testPSM = True
                origG, tri = utils.buildDelaunayGraph(nodes)
                while testPSM:
                    G,tri = utils.buildDelaunayGraph(nodes)
                    d = G.getNodes()
                    psm = 0.0
                    for n in d.iterkeys():
                        node = d[n]
                        for k in node.getNeighbors().iterkeys():
                            other = G.getNode(k)
                            dist = max(utils.lengthBetweenNodes(node,other), 0.000001)
                            try:
                                fac = float(node.width + other.width )/abs(node.x - other.x) #abs(node.x - other.x)#(node.radius + other.radius)/dist
                            except ZeroDivisionError:
                                fac = 1.5
                            try:
                                hfac = float(node.height + other.height)/abs(node.y - other.y) #abs(node.y - other.y)
                            except ZeroDivisionError:
                                hfac = 1.5
                            fac = min(float(fac),hfac)
                            tij = max(float(fac), 1.0)
                            smax = 1.5
                            sij = min(tij, smax)
                            dij = sij*dist
                            if dist >0:
                                newX = node.x + (other.x-node.x) / dist * dij
                                newY = node.y + (other.y-node.y) / dist * dij
                                other.disp_x +=newX - other.x
                                other.disp_y += newY - other.y


                            psm += (1/(dij*dij)) * ((dist-dij)*(dist-dij))
                    arcpy.AddMessage("Current psm %s"%(psm))
                    nodes = []
                    for n in d.iterkeys():
                        node = d[n]
                        node.x += node.disp_x
                        node.y += node.disp_y
                        node.disp_x = 0.0
                        node.disp_y = 0.0
                        node.clearNeighbors()
                        nodes.append(node)
                    if psm >= 0.00000000001:
                        testPSM = True
                    else:
                        testPSM = False
                #check for overlaps amongst everything
                testPSM = True
                while testPSM:
                    G, tri = utils.buildDelaunayGraph(nodes)
                    d = G.getNodes()
                    psm = 0.0
                    for n in d.iterkeys():
                        node = d[n]
                        #nodeExt = utils.extentToPolygon(node)
                        for k in d.iterkeys():
                            if n!=k:
                                other = d[k]
                                if node.geometry and other.geometry:
                                    if node.geometry.overlaps(other.geometry):
                                        G.addUndirectedEdge(node.id, other.id, utils.lengthBetweenNodes(node, other))
                    d = G.getNodes()
                    for n in d.iterkeys():
                        node = d[n]
                        for k in node.getNeighbors().iterkeys():
                            other = G.getNode(k)
                            dist = max(utils.lengthBetweenNodes(node,other), 0.000001)
                            try:
                                fac = float(node.width + other.width )/abs(node.x - other.x) #abs(node.x - other.x)#(node.radius + other.radius)/dist
                            except ZeroDivisionError:
                                fac = 1.5
                            try:
                                hfac = float(node.height + other.height)/abs(node.y - other.y) #abs(node.y - other.y)
                            except ZeroDivisionError:
                                hfac = 1.5
                            fac = min(float(fac),hfac)
                            tij = max(float(fac), 1.0)
                            smax = 1.5
                            sij = min(tij, smax)
                            dij = sij*dist
                            if dist >0:
                                newX = node.x + (other.x-node.x) / dist * dij
                                newY = node.y + (other.y-node.y) / dist * dij
                                other.disp_x +=newX - other.x
                                other.disp_y += newY - other.y
                            psm += (float(1)/(dij*dij)) * ((dist-dij)*(dist-dij))
                    arcpy.AddMessage("Current psm %s"%(psm))
                    nodes = []
                    for n in d.iterkeys():
                        node = d[n]
                        node.x += node.disp_x
                        node.y += node.disp_y
                        node.disp_x = 0.0
                        node.disp_y = 0.0
                        node.clearNeighbors()
                        nodes.append(node)
                    if psm >= 0.00000000001:
                        testPSM = True
                    else:
                        testPSM = False
                ic = arcpy.da.InsertCursor(outFC, ["SHAPE@",idFld.name, valFld.name])

                for n in nodes:
                    if not outline:
                        if n.geometry:
                            n.geometry = utils.moveGeometry(n.geometry, n.o_x-n.x, n.o_y-n.y)
                            if n.geometry:
                                ic.insertRow([n.geometry, n.altid,n.value])
                    else:
                        arr = arcpy.Array()
                        arr.add(arcpy.Point(n.x - n.width, n.y-n.height))
                        arr.add(arcpy.Point(n.x - n.width, n.y+n.height))
                        arr.add(arcpy.Point(n.x + n.width, n.y+n.height))
                        arr.add(arcpy.Point(n.x + n.width, n.y-n.height))
                        arr.add(arcpy.Point(n.x - n.width, n.y-n.height))
                        ic.insertRow([arcpy.Polygon(arr), n.altid,n.value])
                del ic
        except:
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
            arcpy.AddError("While creating a Noncontiguous Cartogram the following errors occured:")
            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)







class relationObject(object):
    def __init__(self, ident, value):
        self.ident = ident
        self.orig_x = 0.0
        self.orig_y = 0.0
        self.value = value
        self.boundaries = {}
        self.radius = 0.0
        self.new_x = 0.0
        self.new_y = 0.0
        self.perimeter = 0.0
        self.buffer = arcpy.Geometry()
        self.disp_x = 0.0
        self.disp_y = 0.0
    def addBoundary(self,ident,length):
        if not self.boundaries.has_key(ident):
            self.boundaries[ident]=length

    def getArcPointGeometry(self):
        return arcpy.PointGeometry(arcpy.Point(self.new_x, self.new_y))


    #-------------------------------------------------------------------------------
# Name:        delauney triangulation
# Purpose:
#
#   Python Port of Paul Bourke original C implementation, and Morten Nielsen's
#   Dot Net port  http://paulbourke.net/papers/triangulate/
# Author:      David S Lamb
#
# Created:     11/1/2012
# Copyright:   (c) dslamb 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
class Delaunay(object):
    @staticmethod
    def triangulatePoints(pnts):
        """pnts - list of nodes"""
        i = 0
        j = 0
        nv = len(pnts)

        if nv < 3: return None

        trimax = 4 * nv;

        xmin = pnts[0].x
        ymin = pnts[0].y
        xmax = xmin
        ymax = ymin

        vertex = {}

        for p in pnts:
            vertex[p.id] = p

            if (vertex[p.id].x < xmin): xmin = vertex[p.id].x
            if (vertex[p.id].x > xmax): xmax = vertex[p.id].x
            if (vertex[p.id].y < ymin): ymin = vertex[p.id].y
            if (vertex[p.id].y > ymax): ymax = vertex[p.id].y


        i = 0

        dx = xmax - xmin
        dy = ymax - ymin
        dmax = 0
        if dx>dy:
            dmax = dx
        else:
            dmax = dy

        xmid = (xmax + xmin)*0.5
        ymid = (ymax + ymin) * 0.5

        #create super triangle
        curID = nv
        vertex[curID] = Node(curID, 0, curID)
        vertex[curID].x = (xmid - 2 * dmax)
        vertex[curID].y = (ymid - dmax)


        curID = nv + 1
        vertex[curID] = Node(curID, 0, curID)
        vertex[curID].x = xmid
        vertex[curID].y = (ymid + 2 * dmax)


        curID = nv + 2
        vertex[curID] = Node(curID, 0, curID)
        vertex[curID].x = (xmid + 2 * dmax)
        vertex[curID].y = (ymid - dmax)




        triangles = []
        triangles.append(DelaunayTriangle(vertex[nv],vertex[nv+1],vertex[nv+2]))
        #add each point to existing mesh
        i=0
        for i in range(nv):
            dEdges = []
            j = 0

            for j in range(len(triangles)-1,-1,-1):
                try:
                    t=triangles[j]
                    if (Delaunay.InCircle(vertex[i], triangles[j].node1, triangles[j].node2, triangles[j].node3)):
                        dEdges.append(DelaunayEdge(t.node1, t.node2))
                        dEdges.append(DelaunayEdge(t.node2, t.node3))
                        dEdges.append(DelaunayEdge(t.node3, t.node1))
                        del triangles[j]
                except:
                    #arcpy.AddMessage( "out of range")
                    continue

            if i >= nv:
                continue


            j = len(dEdges) - 2

            toRemove ={}

            for j in range(len(dEdges)-2,-1,-1):
                for k in range(len(dEdges)-1,j,-1):
                    if j!=k:
                        if dEdges[j]==dEdges[k]:
                            if not toRemove.has_key(j):
                                toRemove[j]=True
                            if not toRemove.has_key(k):
                                toRemove[k]=True

            for j in sorted(toRemove.iterkeys(), reverse=True):
                try:
                    del dEdges[j]
                except Exception as err:
                    print err
                    continue



            j = 0
            for j in range(0,len(dEdges),1):
                if(len(triangles)>=trimax):
                   pass #arcpy.AddMessage( "Exceed maximum Delaunay Edges")
                triangles.append(DelaunayTriangle(dEdges[j].startNode, dEdges[j].endNode, vertex[i]))
            dEdges = None


        i = len(triangles) - 1
        while i >= 0:
            if (triangles[i].node1.id >= nv or triangles[i].node2.id >= nv or triangles[i].node3.id>=nv):
                del triangles[i]
            i-=1
        del vertex[nv]
        del vertex[nv+1]
        del vertex[nv+2]

        return triangles

    @staticmethod
    def InCircle(node, node1, node2, node3):
        Epsilon = .0000000000001

        if (math.fabs(node1.y - node2.y)<Epsilon and math.fabs(node2.y - node3.y)<Epsilon):
            return False
        m1 = 0.0
        m2 = 0.0
        mx1 = 0.0
        mx2 = 0.0
        my1 = 0.0
        my2 = 0.0
        xc = 0.0
        yc = 0.0

        if(math.fabs(node2.y - node1.y)<Epsilon):
            print "if %s - %s" %(node2.y,node1.y)
            m2 = -(node3.x - node2.x) / (node3.y - node2.y)
            mx2 = (node2.x + node3.x)*0.5
            my2 = (node2.y + node3.y)*0.5

            xc = (node2.x + node1.x) * 0.5
            yc = m2 * (xc - mx2) + my2
        elif (math.fabs(node3.y - node2.y)<Epsilon):
            print "el if %s - %s" %(node3.y,node2.y)
            m1 =-(node2.x - node1.x) /(node2.y - node1.y)
            mx1 = (node1.x + node2.x) * 0.5
            my1 = (node1.y + node2.y) * 0.5
            xc = (node3.x + node2.x) * 0.5
            yc = m1 * (xc - mx1) + my1
        else:
            print "else %s - %s" %(node2.y,node1.y)
            m1 = -(node2.x - node1.x) / (node2.y - node1.y)
            m2 = -(node3.x - node2.x) / (node3.y - node2.y)
            mx1 = (node1.x + node2.x) * 0.5
            mx2 = (node2.x + node3.x) * 0.5
            my1 = (node1.y + node2.y) * 0.5
            my2 = (node2.y + node3.y) * 0.5
            xc = (m1 * mx1 - m2 * mx2 + my2 - my1) /(m1-m2)
            yc = m1 * (xc - mx1) + my1

        dx = node2.x - xc
        dy = node2.y - yc
        rsqr = dx * dx + dy * dy
        dx = node.x - xc
        dy = node.y - yc
        drsqr = dx * dx + dy * dy
        if drsqr <= rsqr:
            return True
        return False
    @staticmethod
    def minFloat(x):
        while(x*0.5 != 0):
            x = x*0.5
        return x

class DelaunayTriangle(object):
    def __init__(self, node1, node2, node3):
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.center = Node(0,0,0)
        self.mid0 = Node(1,0,1)
        self.mid1 = Node(2,0,2)
        self.mid2 = Node(3,0,3)
        #self.cirlce = circle()
        self.getCenter()
        self.getEdgeCenter()


    def getCenter(self):
        if self.center:
            self.center.x = (self.node1.x +self.node2.x + self.node3.x)/3
            self.center.y = (self.node1.y +self.node2.y + self.node3.y)/3

    def getEdgeCenter(self):
        self.mid0.x = self.node1.x + (self.node2.x - self.node1.x)/2
        self.mid0.y = self.node1.y + (self.node2.y - self.node1.y)/2

        self.mid1.x = self.node2.x + (self.node3.x - self.node2.x)/2
        self.mid1.y = self.node2.y + (self.node3.y - self.node2.y)/2

        self.mid2.x = self.node3.x + (self.node1.x - self.node3.x)/2
        self.mid2.y = self.node3.y + (self.node1.y - self.node3.y)/2

    def nodePartOfTriangle(self, node):
        if self.node1.id == node.id:
            return True
        if self.node2.id == node.id:
            return True
        if self.node3.id == node.id:
            return True

class DelaunayEdge(object):
    def __init__(self, startNode, endNode):
        self.startNode = startNode
        self.endNode = endNode

    def __eq__(self, dEdge):
        if dEdge.startNode == self.endNode and dEdge.endNode == self.startNode or dEdge.startNode == self.startNode and dEdge.endNode == self.endNode:
                return True
        return False



###Graph classes
class Node(object):
    def __init__(self,id, value, altid):
        self.id = id
        self.altid = altid
        self.neighbors = {}
        self.neighbors_naturalLength = {}
        self.neighbors_springStrength = {}
        self.neighbors_force = {}
        self.degree= 0
        self.x= 0.0
        self.y = 0.0
        self.o_x = 0.0
        self.o_y = 0.0
        self.forceR =[0.0,0.0]
        self.forceA =[0.0,0.0]
        self.velocity =[0.0,0.0]
        self.radius=0.0
        self.value = value
        self.sqrt_value = math.sqrt(self.value)
        self.disp_x = 0.0
        self.disp_y = 0.0
        self.scale_fact = 1.0
        self.geometry = arcpy.Geometry()
        self.width = 0.0
        self.height = 0.0
        self.area = 0.0

    def getNeighbors(self):
        return self.neighbors

    def addNeighbor(self, node, edgeLength):
        if not self.neighbors.has_key(node.id):
            self.neighbors[node.id] = edgeLength
            self.neighbors_naturalLength[node.id] = 0.0
            self.neighbors_springStrength[node.id]  = 0.0
            self.neighbors_force[node.id]=0.0
            self.degree +=1
        else:
            self.neighbors[node.id] = edgeLength

    def deleteNeighbor(self, node):
        if self.neighbors.has_key(node.id):
            del self.neighbors[node.id]
            del self.neighbors_naturalLength[node.id]
            del self.neighbors_springStrength[node.id]
            del self.neighbors_force[node.id]
            self.degree -= 1

    def getLengthToNeighbor(self, neighbor):
        if self.neighbors.has_key(neighbor.id):
            return self.neighbors[neighbor.id]
        return None

    def clearNeighbors(self):
        self.neighbors = {}


class Graph(object):
    def __init__(self):
        self.nodes = {}

    def addNode(self, id, value=0, altid=0):
        if not self.nodes.has_key(id):
            self.nodes[id]=Node(id,value, altid)
        return self.nodes[id]

    def addNodeClass(self, node):
        if not self.nodes.has_key(node.id):
            self.nodes[node.id]=node

    def addUndirectedEdge(self, id1, id2, length):
        if id1 != id2:
            self.nodes[id1].addNeighbor(self.nodes[id2], length)
            self.nodes[id2].addNeighbor(self.nodes[id1], length)
        return

    def addDirectedEdge(self, id1, id2, length):
        if id1 != id2:
            self.nodes[id1].addNeighbor(self.nodes[id2], length)
        return

    def deleteUndirectedEdge(self, id1, id2):
        self.nodes[id1].deleteNeighbor(self.nodes[id2])
        self.nodes[id2].deleteNeighbor(self.nodes[id1])

    def nodeIDList(self):
        listIDs = []
        for k in self.nodes:
            listIDs.append(k)

    def getNodes(self):
        return self.nodes

    def getNode(self, id):
        if self.nodes.has_key(id):
            return self.nodes[id]
    def edgeTupleList(self):
        return

    def getNodeCount(self):
        return len(self.nodes)

    def removeNode(self, nodeID):
        if self.nodes.has_key(nodeID):
            del self.nodes[nodeID]
            for k in self.nodes.iterkeys():
                n = self.nodes[k]
                if n.getNeighbors().has_key(nodeID):
                    del n.getNeighbors()[nodeID]
        return



class utils(object):
    @staticmethod
    def lengthBetweenNodes(node1, node2):
        dx= node1.x - node2.x
        dy= node1.y - node2.y
        return math.sqrt(dx*dx+dy*dy)
    @staticmethod
    def getBearingRadians(node1, node2):
        l = utils.lengthBetweenNodes(node1,node2)
        dx= node1.x - node2.x
        dy = node1.y - node2.y
        return math.atan(dx/dy)
    @staticmethod
    def getDeltaXY(distance, angRadians):
        """Returns a tuple of delta X and Y"""
        return ((distance*math.sin(angRadians)),(distance*math.cos(angRadians)))

    @staticmethod
    def resizeGeometry(geometry, scaleFact, minArea):
        """Scale Geometry by Scale Factor
            geometry = arcpy.Geometry
            scaleFact = floating point scale factor
            minArea = parts that have an area less than this value will not be scaled
            returns a new scaled geometry"""
        typ = string.upper(geometry.type)
        if typ == "POLYGON" or typ == "POLYLINE":
            partsArray = arcpy.Array()
            center = geometry.centroid
            for part in geometry:
                temp = arcpy.Polygon(part)
                pntArray = arcpy.Array()
                for pnt in part:
                    if pnt:
                        dist =(pnt.X - center.X, pnt.Y - center.Y)
                        newPoint = arcpy.Point(dist[0] * scaleFact + center.X, dist[1] * scaleFact + center.Y)
                        pntArray.add(newPoint)
                temp = arcpy.Polygon(pntArray)
                if temp.area >= minArea:
                    partsArray.add(pntArray)
            if partsArray.count > 0:
                if typ == "POLYGON":
                    return arcpy.Polygon(partsArray)
                if type == "POLYLINE":
                    return arcpy.Polyline(partsArray)
        return None

    @staticmethod
    def moveGeometry(geometry, xDisp, yDisp):
        """Shifts polygon and polyline geometry by the x and y displacement
        returns new geometry"""
        typ = string.upper(geometry.type)
        if typ == "POLYGON" or typ == "POLYLINE":
            partsArray = arcpy.Array()
            for part in geometry:
                pntArray = arcpy.Array()
                for pnt in part:
                    newPoint = arcpy.Point(pnt.X - xDisp, pnt.Y - yDisp)
                    pntArray.add(newPoint)
                partsArray.add(pntArray)
            if partsArray.count > 0:
                if typ == "POLYGON":
                    return arcpy.Polygon(partsArray)
                if type == "POLYLINE":
                    return arcpy.Polyline(partsArray)
        return None
    @staticmethod
    def extentToPolygon(n):
        arr = arcpy.Array()
        hw = float(n.geometry.extent.width) / 2
        hh = float(n.geometry.extent.height) / 2
        arr.add(arcpy.Point(n.x - hw, n.y-hw))
        arr.add(arcpy.Point(n.x - hw, n.y+hw))
        arr.add(arcpy.Point(n.x + hw, n.y+hw))
        arr.add(arcpy.Point(n.x + hw, n.y-hw))
        arr.add(arcpy.Point(n.x - hw, n.y-hw))
        return arcpy.Polygon(arr)

    @staticmethod
    def buildPointQuadTree(pnts):
        """List of ArcPy Points"""
        xMax = 0.0
        xMin = math.pow(10,10)
        yMax = 0.0
        yMin = math.pow(10,10)

        for pnt in pnts:
            if pnt.X > xMax:
                xMax = pnt.X
            if pnt.X < xMin:
                xMin = pnt.X
            if pnt.Y > yMax:
                yMax = pnt.Y
            if pnt.Y < yMin:
                yMin = pnt.Y
        lleft = arcpy.Point(xMin, yMin)
        uleft = arcpy.Point(xMin, yMax)
        uright = arcpy.Point(xMax, yMax)
        lright = arcpy.Point(xMax, yMin)

        zerouleft = arcpy.Point(xMin, ((yMax - yMin)*.5 + yMin))
        zerolright = arcpy.Point(((xMax - xMin)*.5 + xMin),yMin)
        zerouright = arcpy.Point(((xMax - xMin)*.5 + xMin), ((yMax - yMin)*.5 + yMin))

        twouright = arcpy.Point(((xMax - xMin)*.5 + xMin), yMax)
        oneuright = arcpy.Point(xMax, ((yMax - yMin)*.5 + yMin))

        G = Graph()
        n=G.addNode("root",0,-1)
        n.geometry = utils.build4Polygon(arcpy.Point(xMin, yMin), arcpy.Point(xMin, yMax), arcpy.Point(xMax, yMax), arcpy.Point(xMax, yMin))

        n = G.addNode('0',0,0)
        n.geometry = utils.build4Polygon(lleft, zerouleft, zerouright, zerolright)
        G.addDirectedEdge("root", '0', 0)
        utils.recurseSplitQuads(G, n, pnts)
        n = G.addNode('1',0,1)
        n.geometry = utils.build4Polygon(zerolright, zerouright, oneuright, lright)
        G.addDirectedEdge("root", '1', 0)
        utils.recurseSplitQuads(G, n, pnts)
        n = G.addNode('2',0,1)
        n.geometry = utils.build4Polygon(zerouleft, uleft, twouright, zerouright)
        G.addDirectedEdge("root", '2', 0)
        utils.recurseSplitQuads(G, n, pnts)
        n = G.addNode('3',0,1)
        n.geometry = utils.build4Polygon(zerouright, twouright, uright, oneuright)
        G.addDirectedEdge("root", '3', 0)
        utils.recurseSplitQuads(G, n, pnts)

        return G



    @staticmethod
    def recurseSplitQuads(G, node, pnts):
        containerCount = 0
        duplicate = {}
        for pnt in pnts:
            t = "%s|%s"%(pnt.X, pnt.Y)
            if not duplicate.has_key(t):
                if node.geometry.contains(pnt):
                    containerCount += 1
                    duplicate[t] = True
        if containerCount > 1:
            xMax = node.geometry.extent.XMax
            xMin = node.geometry.extent.XMin
            yMax = node.geometry.extent.YMax
            yMin = node.geometry.extent.YMin
            lleft = arcpy.Point(xMin, yMin)
            uleft = arcpy.Point(xMin, yMax)
            uright = arcpy.Point(xMax, yMax)
            lright = arcpy.Point(xMax, yMin)

            zerouleft = arcpy.Point(xMin, ((yMax - yMin)*.5 + yMin))
            zerolright = arcpy.Point(((xMax - xMin)*.5 + xMin),yMin)
            zerouright = arcpy.Point(((xMax - xMin)*.5 + xMin), ((yMax - yMin)*.5 + yMin))

            twouright = arcpy.Point(((xMax - xMin)*.5 + xMin), yMax)
            oneuright = arcpy.Point(xMax, ((yMax - yMin)*.5 + yMin))

            arcpy.AddMessage("\n"+node.id+"\n")

            n = G.addNode(node.id + '0',0,0)
            n.geometry = utils.build4Polygon(lleft, zerouleft, zerouright, zerolright)
            arcpy.AddMessage(n.id)
            G.addDirectedEdge(node.id, n.id, 0)
            utils.recurseSplitQuads(G,n, pnts)

            n = G.addNode(node.id + '1',0,1)
            n.geometry = utils.build4Polygon(zerolright, zerouright, oneuright, lright)
            arcpy.AddMessage(n.id)
            G.addDirectedEdge(node.id, n.id, 0)
            utils.recurseSplitQuads(G, n, pnts)

            n = G.addNode(node.id + '2',0,1)
            n.geometry = utils.build4Polygon(zerouleft, uleft, twouright, zerouright)
            arcpy.AddMessage(n.id)
            G.addDirectedEdge(node.id, n.id, 0)
            utils.recurseSplitQuads(G, n, pnts)

            n = G.addNode(node.id + '3',0,1)
            arcpy.AddMessage(n.id)
            n.geometry = utils.build4Polygon(zerouright, twouright, uright, oneuright)
            G.addDirectedEdge(node.id, n.id, 0)
            utils.recurseSplitQuads(G, n, pnts)
        else:
            return True
        return True

    @staticmethod
    def build4Polygon(pnt1, pnt2, pnt3, pnt4):
        try:
            arr = arcpy.Array()
            arr.add(pnt1)
            arr.add(pnt2)
            arr.add(pnt3)
            arr.add(pnt4)
            arr.add(pnt1)
            return arcpy.Polygon(arr)
        except:
            return none


    @staticmethod
    def minimalSpanningTree(G, startNode):
        """From Python Algorithms by Magnus Lie Hetland.  Listing 7-5 Prim's Algorithm page 167"""
        previous = {}
        queue = [(0,None,startNode.id)]
        while queue:
            _, p, n = heappop(queue)
            if n in previous: continue
            previous[n] = p
            node = G.getNode(n)
            print node

            if node:
                print node.id
                neighbors = node.getNeighbors()
                for k in neighbors.iterkeys():
                    length = node.getNeighbors()[k]
                    heappush(queue,(length, n, k.id))
        spanningG = Graph()
        for id in G.getNodes().iterkeys():
            spanningG.addNode(id)
        for k in previous.iterkeys():
            p = previous[k]
            if p:
                spanningG.addDirectedEdge(p, k, utils.lengthBetweenNodes(spanningG.getNode(p), spanningG.getNode(k)))

        return previous, spanningG

    @staticmethod
    def Dijkstra(G,start,end=None):
        D = {}   # dictionary of final distances
        P = {}   # dictionary of previous nodes
        Q = [(0,None,start)]  # heapq for tuple(distance, predecessors, and node)...start with first node, no predecessors
        while Q:
            dist,pred,v = heappop(Q)  #pop priority node...distance sets priority
            if v in D:
                continue  # Check that a node distance hasn't been added yet
            D[v] = dist  #node distance
            P[v] = pred
            cNeighbors = G.getNode(v.id).getNeighbors()
            for w in cNeighbors:  #loop through neighbors
                heappush(Q, (D[v] + cNeighbors[w], v, G.getNode(w.id)))  #push priority distance
        return (D,P)  #return distances and predecessors.  You use the predecessors to get path.

    @staticmethod
    def shortestPath(G,start,end):
    	D,P = utils.Dijkstra(G,start,end)  #calculate the shortest path, get distance and nodes
    	Path = []
    	while 1:
    		Path.append(end.id)
    		if end == start: break
    		end = P[end]
    	Path.reverse()
    	return Path

    @staticmethod
    def buildDelaunayGraph(nodes):
        G = Graph()
        triangles = Delaunay.triangulatePoints(nodes)
        for n in nodes:
            G.addNodeClass(n)
        for e in triangles:
            G.addNodeClass(e.node1)
            G.addNodeClass(e.node2)
            G.addNodeClass(e.node3)
            G.addUndirectedEdge(e.node1.id, e.node2.id, utils.lengthBetweenNodes(e.node1, e.node2))
            G.addUndirectedEdge(e.node2.id, e.node3.id, utils.lengthBetweenNodes(e.node2, e.node3))
            G.addUndirectedEdge(e.node3.id, e.node1.id, utils.lengthBetweenNodes(e.node3, e.node1))
        return G,triangles

    @staticmethod
    def getBisector(node1, node2, length):
        """Perpendicular bisector between two nodes, returns tuple of arcpy points"""
        m = float(node1.y - node2.y)/float(node1.x - node2.x)
        invM = float(-1)/m
        midX = float(node1.x + node2.x)/2
        midY = float(node1.y + node2.y)/2
        newNode = arcpy.Point()
        newNode2 = arcpy.Point()
        newNode.X = midX + (float(length) / math.sqrt(1 + invM*invM))
        newNode.Y = midY + (invM * length / math.sqrt(1+invM*invM))
        newNode2.X = midX - (float(length) / math.sqrt(1 + invM*invM))
        newNode2.Y = midY - (invM * length / math.sqrt(1+invM*invM))
        return (newNode,newNode2)

    @staticmethod
    def pointIntersection(node11, node12, node21, node22):
        """arcpy points  returns arcpy point, or None"""
        m1 = (node11.Y - node12.Y)/(node11.X - node12.X)
        m2 = (node21.Y - node22.Y)/(node21.X - node22.X)
        if (m1 != m2):
            if m1 and m2:
                b1 = node11.Y -(m1 * node11.X)
                b2 = node21.Y - (m2 * node21.X)
                x = (b2 - b1) / (m1 - m2)
                y = (m1 * x) + b1
                return arcpy.Point(x,y)
        else:
            return None

    @staticmethod
    def circleFromTriangle(dt1):
        """delaunay triangles, returns circle class"""
        ax = dt1.node1.x
        ay = dt1.node1.y
        bx = dt1.node2.x
        by = dt1.node2.y
        cx = dt1.node3.x
        cy = dt1.node3.y
        D = 2.0 * (ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
        a = utils.lengthBetweenNodes(dt1.node1, dt1.node2)
        b = utils.lengthBetweenNodes(dt1.node2, dt1.node3)
        c = utils.lengthBetweenNodes(dt1.node3, dt1.node1)
        circ = circle()
        circ.centerX = float((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by)) / D
        circ.centerY = float((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax)) / D
        circ.radius = float(a*b*c)/math.sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c))
        return circ

    @staticmethod
    def voronoiPolygons(pnts,ext):
        """list of nodes (node class), arcpy extent object, returns list of pnts that create a polygon"""
        triangles = Delaunay.triangulatePoints(pnts)
        polygons = {}
        for n in pnts:
            basePoly = []
            for t in triangles:
                if t.nodePartOfTriangle(n):
                    c = utils.circleFromTriangle(t)
                    pnt = arcpy.Point(c.centerX, c.centerY)
                    #if ext.contains(pnt):
                    basePoly.append(pnt)
            angles = []
            if len(basePoly)<2:
                continue
            for j in basePoly:
                angles.append(math.atan2(j.Y-n.y, j.X - n.x) + math.pi * 2)

            for j in range(0,len(basePoly)-1):
                for k in range(0, len(basePoly)-1-j):
                    if (angles[k] > angles[k+1]):
                        p1 = basePoly[k]
                        basePoly[k] = basePoly[k+1]
                        basePoly[k+1] = p1
                        tangle = angles[k]
                        angles[k] = angles[k+1]
                        angles[k+1] = tangle

            pArray = arcpy.Array()
            for j in basePoly:
                pArray.add(j)
            pgon = arcpy.Polygon(pArray)
            clip = pgon.clip(ext)
            polygons[n.id] = clip
        return polygons

    @staticmethod
    def signum(d):
        if d == 0.0:
            return 0.0
        if d > 0.0:
            return 1.0
        if d < 0.0:
            return -1.0