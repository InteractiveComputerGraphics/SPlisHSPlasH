import maya.cmds as cmds
import maya.api.OpenMaya as apiOpenMaya
import maya.OpenMayaMPx as OpenMayaMPx
import json
import os
import math
import sys
import re
import struct
from collections import OrderedDict
from copy import deepcopy
import maya.OpenMaya as OpenMaya
import maya.OpenMayaMPx as OpenMayaMPx


	
######################################################
# PluginFunctions
######################################################
class PluginFunctions():

	######################################################
	# getAllNodesOfType
	######################################################	
	@staticmethod	
	def getAllNodesOfType(typeId):
		list = cmds.ls( type='transform', long=True ) 
		result = []
		for node in list:
			# find type attribute
			sphAttr = cmds.listAttr(node, string="SPH_Type")
			if sphAttr != None:
				sphtype = cmds.getAttr(node + ".SPH_Type")
				if typeId == sphtype:
					result.append(node)
		return result
		
	######################################################
	# getShape
	######################################################	
	@staticmethod
	def getShape(nodeName):
		return cmds.listRelatives(nodeName, shapes=True, type="shape")
		
	######################################################
	# get quaternion of a transform node
	######################################################
	@staticmethod
	def getQuaternion(node):
		sel_list = apiOpenMaya.MSelectionList()
		sel_list.add(node)
		obj = sel_list.getDependNode(0)
		xform = apiOpenMaya.MFnTransform(obj)
		quat = xform.rotation(asQuaternion=True)
		quat.normalizeIt()
		return quat
		
	######################################################
	# get axis,angle of a transform node
	######################################################
	@staticmethod
	def getAxisAngle(node):
		sel_list = apiOpenMaya.MSelectionList()
		sel_list.add(node)
		obj = sel_list.getDependNode(0)
		xform = apiOpenMaya.MFnTransform(obj)
		quat = xform.rotation(asQuaternion=True)
		quat.normalizeIt()
		aa = quat.asAxisAngle()
		return ([aa[0][0], aa[0][1], aa[0][2]], aa[1])
	
	@staticmethod
	def createFloatAttr(longName, shortName, defaultValue, softMin, softMax, minValue=0, maxValue=1000000):
		nAttr = OpenMaya.MFnNumericAttribute()
		newAttr = nAttr.create( longName, shortName, OpenMaya.MFnNumericData.kFloat, defaultValue )
		nAttr.setStorable(1)
		nAttr.setMin(minValue)
		nAttr.setMax(maxValue)
		nAttr.setSoftMin(softMin)
		nAttr.setSoftMax(softMax)	
		return newAttr
		
	@staticmethod
	def createIntAttr(longName, shortName, defaultValue, softMin, softMax, minValue=0, maxValue=1000000):
		nAttr = OpenMaya.MFnNumericAttribute()
		newAttr = nAttr.create( longName, shortName, OpenMaya.MFnNumericData.kInt, defaultValue )
		nAttr.setStorable(1)
		nAttr.setMin(minValue)
		nAttr.setMax(maxValue)
		nAttr.setSoftMin(softMin)
		nAttr.setSoftMax(softMax)	
		return newAttr
		
	@staticmethod
	def createBoolAttr(longName, shortName, defaultValue):
		nAttr = OpenMaya.MFnNumericAttribute()
		newAttr = nAttr.create( longName, shortName, OpenMaya.MFnNumericData.kBoolean, defaultValue )
		nAttr.setStorable(1)	
		return newAttr
		
	@staticmethod
	def createVec3Attr(longName, shortName, defaultValue):
		nAttr = OpenMaya.MFnNumericAttribute()
		newAttr = nAttr.create( longName, shortName, OpenMaya.MFnNumericData.k3Float )
		nAttr.setDefault(defaultValue[0], defaultValue[1], defaultValue[2])
		nAttr.setStorable(1)	
		return newAttr
		
	@staticmethod
	def createEnumAttr(longName, shortName, defaultValue, enumList):
		eAttr = OpenMaya.MFnEnumAttribute()
		newAttr = eAttr.create( longName, shortName, defaultValue)
		i=0
		for item in enumList:
			eAttr.addField(item, i)
			i+=1
		eAttr.setStorable(1)	
		return newAttr
		
	######################################################
	# createBoolParam
	######################################################	
	@staticmethod	
	def createBoolParam(name, label, description, defaultValue):
		param = {
				"type": "bool",
				"name": name,
				"label": label,
				"description": description,
				"default": defaultValue,
				"value": defaultValue,
				"ctrlId": None
			}
		return param

	######################################################
	# createFloatParam
	######################################################	
	@staticmethod		
	def createFloatParam(name, label, description, defaultValue, minValue, maxValue, fieldMin=0, fieldMax=1000000):
		param = {
				"type": "float",
				"name": name,
				"label": label,
				"description": description,
				"default": defaultValue,
				"value": defaultValue,
				"min": minValue,
				"max": maxValue,
				"fieldMin": fieldMin,
				"fieldMax": fieldMax,
				"ctrlId": None
			}
		return param

	######################################################
	# createVec3Param
	######################################################	
	@staticmethod		
	def createVec3Param(name, label, description, defaultValue):
		param = {
				"type": "vec3",
				"name": name,
				"label": label,
				"description": description,
				"default": defaultValue,
				"value": defaultValue,
				"ctrlId": None
			}
		return param	
		
	######################################################
	# createIntParam
	######################################################	
	@staticmethod		
	def createIntParam(name, label, description, defaultValue, minValue, maxValue, fieldMin=0, fieldMax=1000000):
		param = {
				"type": "int",
				"name": name,
				"label": label,
				"description": description,
				"default": defaultValue,
				"value": defaultValue,
				"min": minValue,
				"max": maxValue,
				"fieldMin": fieldMin,
				"fieldMax": fieldMax,
				"ctrlId": None
			}
		return param

	######################################################
	# createStringParam
	######################################################	
	@staticmethod		
	def createStringParam(name, label, description, defaultValue):
		param = {
				"type": "string",
				"name": name,
				"label": label,
				"description": description,
				"default": defaultValue,
				"value": defaultValue,
				"ctrlId": None
			}
		return param	
		
	######################################################
	# createEnumParam
	######################################################		
	@staticmethod		
	def createEnumParam(name, label, description, defaultValue, enumList):
		param = {
				"type": "enum",
				"name": name,
				"label": label,
				"description": description,
				"default": defaultValue,
				"value": defaultValue,
				"enumList": enumList,
				"ctrlId": None
			}
		return param
		
	######################################################
	# getSelectedTransforms
	# get all selected transform nodes recursively
	######################################################	
	@staticmethod
	def getSelectedTransforms():
		list = cmds.ls( selection=True, type='transform', long=True )    
		transformNodes = []
		for item in list:
			transformNodes.append(item)
			children = cmds.listRelatives(item, ad=True, type="transform")
			if children == None:
				continue
			for child in children:
				transformNodes.append(child)
		return transformNodes
			
	
######################################################
# createCircularEmitter
######################################################		
class createCircularEmitterCmd(OpenMayaMPx.MPxCommand):
	s_name = "createCircularEmitter"

	def __init__(self):
		OpenMayaMPx.MPxCommand.__init__(self)

	@staticmethod
	def creator():
		return createCircularEmitterCmd()

	def doIt(self, args):
		self.redoIt()

	def redoIt(self):	
		self.cyl = cmds.polyCylinder(name="CircularEmitter", r=1, h=0.2, sx=20, sy=1, sz=1, ax=[1,0,0], rcp=0, cuv=3, ch=1)
		cmds.delete(ch=True)
		node = self.cyl[0]
		cmds.delete(node + ".f[40:59]") 
		
		cmds.scale(0.5, 0.5, 0.5, self.cyl[0])
		
		# set type 
		cmds.addAttr(node, longName="SPH_Type", niceName="type",dt="string", hidden=True)
		cmds.setAttr((node + '.SPH_Type'), "CircularEmitter", type="string")

		# velocity
		cmds.addAttr(node, longName="SPH_velocity", niceName="velocity", at="float");
		cmds.setAttr((node + '.SPH_velocity'), 1.0)
		
		# start time
		cmds.addAttr(node, longName="SPH_startTime", niceName="start time", at="float");
		cmds.setAttr((node + '.SPH_startTime'), 0.0)
		
		# velocity
		cmds.addAttr(node, longName="SPH_endTime", niceName="end time", at="float");
		cmds.setAttr((node + '.SPH_endTime'), 100000.0)
		
		# fluid id
		cmds.addAttr(node, longName="SPH_fluidId", niceName="Fluid id", dt="string")
		cmds.setAttr((node + '.SPH_fluidId'), "Fluid", type="string")
			
	def undoIt(self):
		pass

	def isUndoable(self):
		return True
		
		
######################################################
# RectangularEmitter
######################################################		
class createRectangularEmitterCmd(OpenMayaMPx.MPxCommand):
	s_name = "createRectangularEmitter"

	def __init__(self):
		OpenMayaMPx.MPxCommand.__init__(self)

	@staticmethod
	def creator():
		return createRectangularEmitterCmd()

	def doIt(self, args):
		self.redoIt()


	def redoIt(self):	
		self.cube = cmds.polyCube(name="RectangularEmitter", w=0.2, h=1, d=1, sx=1, sy=1, sz=1, ch=1)
		cmds.delete(ch=True)
		node = self.cube[0]
		cmds.delete(node + ".f[4]") 
		
		# set type 
		cmds.addAttr(node, longName="SPH_Type", niceName="type",dt="string", hidden=True)
		cmds.setAttr((node + '.SPH_Type'), "RectangularEmitter", type="string")

		# velocity
		cmds.addAttr(node, longName="SPH_velocity", niceName="velocity", at="float");
		cmds.setAttr((node + '.SPH_velocity'), 1.0)
		
		# start time
		cmds.addAttr(node, longName="SPH_startTime", niceName="start time", at="float");
		cmds.setAttr((node + '.SPH_startTime'), 0.0)
		
		# velocity
		cmds.addAttr(node, longName="SPH_endTime", niceName="end time", at="float");
		cmds.setAttr((node + '.SPH_endTime'), 100000.0)
		
		# fluid id
		cmds.addAttr(node, longName="SPH_fluidId", niceName="Fluid id", dt="string")
		cmds.setAttr((node + '.SPH_fluidId'), "Fluid", type="string")
		
	def undoIt(self):
		pass

	def isUndoable(self):
		return True	
		
######################################################
# AnimationField
######################################################		
class createAnimationFieldCmd(OpenMayaMPx.MPxCommand):
	s_name = "createAnimationField"
	s_shortTypeFlag = '-s'
	s_longTypeFlag = '-shape'

	def __init__(self):
		OpenMayaMPx.MPxCommand.__init__(self)
		
	@staticmethod
	def syntaxCreator():
		syntax = OpenMaya.MSyntax()
		
		syntax.addFlag( createAnimationFieldCmd.s_shortTypeFlag, createAnimationFieldCmd.s_longTypeFlag, OpenMaya.MSyntax.kLong )

		return syntax		

	@staticmethod
	def creator():
		return createAnimationFieldCmd()

	def doIt(self, args):
		argData = OpenMaya.MArgParser( self.syntax(), args )
		self.shapeType = 0
		if argData.isFlagSet( createAnimationFieldCmd.s_shortTypeFlag ):
			self.shapeType = argData.flagArgumentInt(createAnimationFieldCmd.s_shortTypeFlag, 0)
		self.redoIt()


	def redoIt(self):	
		poly = ""
		if self.shapeType == 1:
			poly = cmds.polySphere(name="AnimationField", r=1, sx=20, sy=20, ax=[0,1,0], cuv=2, ch=1)
			cmds.expression(s=poly[0] + ".scaleY=" + poly[0] + ".scaleZ=" + poly[0] + ".scaleX;", o=poly[0])
		elif self.shapeType == 2:
			poly = cmds.polyCylinder(name="AnimationField", r=1, h=1, sx=20, sy=1, ax=[1,0,0], cuv=3, rcp=0, ch=1)
			cmds.expression(s=poly[0] + ".scaleZ=" + poly[0] + ".scaleY;", o=poly[0])
		else:
			poly = cmds.polyCube(name="AnimationField", w=1, h=1, d=1, sx=1, sy=1, sz=1, ch=1)
		cmds.delete(ch=True)
		node = poly[0]
		
		# set type 
		cmds.addAttr(node, longName="SPH_shapeType", niceName="shape type", at="long", hidden=True)
		cmds.setAttr((node + '.SPH_shapeType'), self.shapeType)
		
		# set type 
		cmds.addAttr(node, longName="SPH_Type", niceName="type",dt="string", hidden=True)
		cmds.setAttr((node + '.SPH_Type'), "AnimationField", type="string")
		
		# set particle field 
		cmds.addAttr(node, longName="SPH_particleField", niceName="paricle field",dt="string")
		cmds.setAttr((node + '.SPH_particleField'), "velocity", type="string")

		# set expression
		cmds.addAttr(node, longName="SPH_expressionX", niceName="expression - x",dt="string")
		cmds.setAttr((node + '.SPH_expressionX'), "", type="string")
		
		cmds.addAttr(node, longName="SPH_expressionY", niceName="expression - y",dt="string")
		cmds.setAttr((node + '.SPH_expressionY'), "", type="string")
		
		cmds.addAttr(node, longName="SPH_expressionZ", niceName="expression - z",dt="string")
		cmds.setAttr((node + '.SPH_expressionZ'), "", type="string")
		
	def undoIt(self):
		pass

	def isUndoable(self):
		return True			


######################################################
# convertToFluid
#
# Converts a list of transform nodes to fluid models.
# Only nodes with a shape are converted.
######################################################	
class convertToFluidCmd(OpenMayaMPx.MPxCommand):
	s_name = "convertToFluid"

	def __init__(self):
		OpenMayaMPx.MPxCommand.__init__(self)

	@staticmethod
	def creator():
		return convertToFluidCmd()

	def doIt(self, args):
		self.redoIt()


	def redoIt(self):	
		nodes = PluginFunctions.getSelectedTransforms()
		self.convertToFluid(nodes)
	
	def convertToFluid(self, nodes):
		for node in nodes:
			shapeNode = PluginFunctions.getShape(node)
			if shapeNode != None:
				# remove all existing attributes
				sphAttr = cmds.listAttr(node, string="SPH_*")
				if sphAttr != None:
					for attr in sphAttr:
						cmds.deleteAttr(node, at=attr)
				
				# set type 
				cmds.addAttr(node, longName="SPH_Type", niceName="type",dt="string", hidden=True)
				cmds.setAttr((node + '.SPH_Type'), "Fluid", type="string")

				# fluid id
				cmds.addAttr(node, longName="SPH_fluidId", niceName="Fluid id", dt="string")
				cmds.setAttr((node + '.SPH_fluidId'), "Fluid", type="string")
		
				# initialVelocity
				cmds.addAttr(node, longName="SPH_initialVelocity", niceName="Initial velocity", at="float3")
				cmds.addAttr(node, longName='SPHinitialVelocityX', at='float', parent='SPH_initialVelocity')
				cmds.addAttr(node, longName='SPHinitialVelocityY', at='float', parent='SPH_initialVelocity')
				cmds.addAttr(node, longName='SPHinitialVelocityZ', at='float', parent='SPH_initialVelocity')
				cmds.setAttr((node + '.SPH_initialVelocity'), 0, 0, 0, type="float3")
				
				# resolutionSDF
				cmds.addAttr(node, longName="SPH_resolutionSDF", niceName="SDF resolution", at="long3")
				cmds.addAttr(node, longName='SPHresolutionSDFX', at='long', parent='SPH_resolutionSDF')
				cmds.addAttr(node, longName='SPHresolutionSDFY', at='long', parent='SPH_resolutionSDF')
				cmds.addAttr(node, longName='SPHresolutionSDFZ', at='long', parent='SPH_resolutionSDF')
				cmds.setAttr((node + '.SPH_resolutionSDF'), 20, 20, 20, type="long3")
				
				# invert
				cmds.addAttr(node, longName="SPH_invert", niceName="invert", at="bool");
				cmds.setAttr((node + '.SPH_invert'), False)
				
				# denseMode
				cmds.addAttr(node, longName="SPH_denseMode", niceName="Dense mode", at="enum", enumName="Regular=0:Almost dense=1:Dense=2")
				cmds.setAttr((node + '.SPH_denseMode'), 0)
				
				
######################################################
# convertToRigidBody
#
# Converts a list of transform nodes to rigid bodies.
# Only nodes with a shape are converted.
######################################################		
class convertToRigidBodiesCmd(OpenMayaMPx.MPxCommand):
	s_name = "convertToRigidBodies"

	def __init__(self):
		OpenMayaMPx.MPxCommand.__init__(self)

	@staticmethod
	def creator():
		return convertToRigidBodiesCmd()

	def doIt(self, args):
		self.redoIt()


	def redoIt(self):	
		nodes = PluginFunctions.getSelectedTransforms()
		self.convertToRigidBody(nodes)
	
	def convertToRigidBody(self, nodes):
		for node in nodes:
			shapeNode = PluginFunctions.getShape(node)
			if shapeNode != None:
				# remove all existing attributes
				sphAttr = cmds.listAttr(node, string="SPH_*")
				if sphAttr != None:
					for attr in sphAttr:
						cmds.deleteAttr(node, at=attr)
				
				# set type 
				cmds.addAttr(node, longName="SPH_Type", niceName="type",dt="string", hidden=True)
				cmds.setAttr((node + '.SPH_Type'), "RigidBody", type="string")

				# is dynamic
				cmds.addAttr(node, longName="SPH_isDynamic", niceName="dynamic", at="bool");
				cmds.setAttr((node + '.SPH_isDynamic'), False)
				
				# is wall
				cmds.addAttr(node, longName="SPH_isWall", niceName="wall", at="bool");
				cmds.setAttr((node + '.SPH_isWall'), False)
				
				# is dynamic
				cmds.addAttr(node, longName="SPH_color", niceName="color", usedAsColor= True, at="float3")
				
				# map thickness
				cmds.addAttr(node, longName="SPH_density", niceName="density", at="float");
				cmds.setAttr((node + '.SPH_density'), 1000.0)
				
				# map invert
				cmds.addAttr(node, longName="SPH_mapInvert", niceName="invert map", at="bool");
				cmds.setAttr((node + '.SPH_mapInvert'), False)
				
				# map thickness
				cmds.addAttr(node, longName="SPH_mapThickness", niceName="map thickness", at="float");
				cmds.setAttr((node + '.SPH_mapThickness'), 0.0)
				
				# map resolution
				cmds.addAttr(node, longName="SPH_mapResolution", niceName="map resolution", at="long3")
				cmds.addAttr(node, longName='SPHmapResolutionX', at='long', parent='SPH_mapResolution')
				cmds.addAttr(node, longName='SPHmapResolutionY', at='long', parent='SPH_mapResolution')
				cmds.addAttr(node, longName='SPHmapResolutionZ', at='long', parent='SPH_mapResolution')
				cmds.setAttr((node + '.SPH_mapResolution'), 20, 20, 20, type="long3")
				
						
				cmds.addAttr(node, longName='SPHcolorR', at='float', parent='SPH_color')
				cmds.addAttr(node, longName='SPHcolorG', at='float', parent='SPH_color')
				cmds.addAttr(node, longName='SPHcolorB', at='float', parent='SPH_color')
				
				cmds.setAttr((node + '.SPH_color'), 0.2, 0.2, 0.2, type='float3')

					

######################################################
# saveModel
######################################################					
class saveModelCmd(OpenMayaMPx.MPxCommand):
	s_name = "saveModel"

	def __init__(self):
		OpenMayaMPx.MPxCommand.__init__(self)

	@staticmethod
	def creator():
		return saveModelCmd()

	def doIt(self, args):
		self.redoIt()


	def redoIt(self):	
		sphConfigList = cmds.ls( type='SPHConfigurationNode', long=True ) 
		if len(sphConfigList) == 0:
			cmds.warning("Not saved since no SPH configuration node was found.")
			return 
	
		if not cmds.pluginInfo("objExport", query=True, loaded=True):
			cmds.loadPlugin("objExport")

		fileName = cmds.fileDialog2(ff="*.json", fm=0, dir="")
		scenePath = os.path.dirname(fileName[0])
		
		scene=self.generateScene(scenePath)
		if scene == None:
			return
		
		f = open(fileName[0], 'w')
		json_str = json.dumps(scene, sort_keys=True,indent=4, separators=(',', ': '))
		f.write(json_str)
		f.close()
		
	def isUndoable(self):
		return False
		
	######################################################
	# openFluidIdDialog
	######################################################					
	def openFluidIdDialog(self):
		sphConfigList = cmds.ls( type='SPHConfigurationNode', long=False )
		
		cmds.columnLayout( adjustableColumn=True, columnOffset=["both", 10], rowSpacing=10, columnAlign="center" )
			
		cmds.textScrollList("SPHFluidIdList", numberOfRows=8, allowMultiSelection=False,
				append=sphConfigList,
				selectItem=sphConfigList[0], showIndexedItem=1)
		cmds.rowLayout(numberOfColumns=2)
		cmds.button("Ok", c='cmds.layoutDialog( dismiss="Ok " + cmds.textScrollList("SPHFluidIdList",q=True,selectItem=True)[0]  )' )
		cmds.button("Cancel", c='cmds.layoutDialog( dismiss="Cancel" )')
		
	######################################################
	# generate scene
	######################################################
	def generateScene(self, scenePath):
		scene = OrderedDict()
		scene['FluidModels'] = []
		scene['RigidBodies'] = []
		scene['Emitters'] = []
		scene['AnimationFields'] = []
		scene['Configuration'] = OrderedDict()
			
		sphConfigList = cmds.ls( type='SPHConfigurationNode', long=True ) 
		sphConfig = ""
		if len(sphConfigList) == 0:
			cmds.warning("Not saved since no SPH configuration node was found.")
			return None
		elif len(sphConfigList) > 1:
			sphConfig = sphConfigList[0]
			res = cmds.layoutDialog(ui=self.openFluidIdDialog)
			if res == "Cancel":
				return None
			else:
				sphConfig = res[3:]
		else:
			sphConfig = sphConfigList[0]
				
			#cmds.warning("More than one SPH configuration node was found using " + sphConfigList[0] + ".")
		
		attributes = cmds.listAttr(sphConfig, string="SPH_*", sn=False)
			
		for attr in attributes:
			if cmds.getAttr(sphConfig + "." + attr, type=True) == "float3":
				value = cmds.getAttr(sphConfig + "." + attr)[0]
			else:
				value = cmds.getAttr(sphConfig + "." + attr)
						
			# avoid to write child attributes
			parent = cmds.attributeQuery( attr, node=sphConfig, listParent=True )
			if parent == None:
				scene["Configuration"][attr[4:]] = value
			
		fluidConfigList = cmds.ls( type='SPHFluidConfigurationNode', long=False ) 
		if len(fluidConfigList) == 0:
			cmds.warning("Not saved since no fluid material node was found.")
			return 		

		for fluid in fluidConfigList:
			attributes = cmds.listAttr(fluid, string="SPH_*", sn=False)
			scene[fluid] = OrderedDict()
			for attr in attributes:
				if cmds.getAttr(fluid + "." + attr, type=True) == "float3":
					value = cmds.getAttr(fluid + "." + attr)[0]
				else:
					value = cmds.getAttr(fluid + "." + attr)
				scene[fluid][attr[4:]] = value
					
		rigidBodies = PluginFunctions.getAllNodesOfType("RigidBody")
		for rb in rigidBodies:
			self.addRigidBody(scene, rb, scenePath)
			
		fluidModels = PluginFunctions.getAllNodesOfType("Fluid")
		for fluid in fluidModels:
			self.addFluid(scene, fluid, scenePath)
			
		emitters = PluginFunctions.getAllNodesOfType("RectangularEmitter")
		for emitter in emitters:
			self.addRectangularEmitter(sphConfig, scene, emitter, scenePath)
		
		emitters = PluginFunctions.getAllNodesOfType("CircularEmitter")
		for emitter in emitters:
			self.addCircularEmitter(sphConfig, scene, emitter, scenePath)
			
		animFields = PluginFunctions.getAllNodesOfType("AnimationField")
		for animField in animFields:
			self.addAnimationField(sphConfig, scene, animField, scenePath)
			
		return scene	
		
	######################################################
	# getCurrentParticleRadius
	######################################################					
	def getCurrentParticleRadius(self, sphConfig):
		return cmds.getAttr(sphConfig + ".particleRadius")
		
	######################################################
	# add rigid bodies
	######################################################
	def addRigidBody(self, scene, rbNode, scenePath):
		
		# export geometry
		cmds.select([rbNode], replace=True)
		polyTri = cmds.polyTriangulate()
		name = cmds.ls( selection=True, type='transform', long=False )[0] 
		fileName = os.path.join(scenePath, "rb_" + name + ".obj")
		cmds.file(fileName, force=True, options="groups=0;ptgroups=0;materials=0;smoothing=0;normals=0", pr=True, exportSelected=True, type="OBJexport")
		cmds.delete(polyTri)

		isDynamic = cmds.getAttr(rbNode + ".SPH_isDynamic")
		isWall = cmds.getAttr(rbNode + ".SPH_isWall")
		density = cmds.getAttr(rbNode + ".SPH_density")
		mapInvert = cmds.getAttr(rbNode + ".SPH_mapInvert")
		mapThickness = cmds.getAttr(rbNode + ".SPH_mapThickness")
		mapResolution = cmds.getAttr(rbNode + ".SPH_mapResolution")[0]
		color = cmds.getAttr(rbNode + ".SPH_color")[0]
		color = color + (1.0,)
		rb = {  'geometryFile': "rb_" + name + ".obj",
				'translation': [0,0,0],
				'rotationAxis': [1,0,0],
				'rotationAngle': 0.0,
				'scale': [1,1,1],
				'color': color,
				'isDynamic': isDynamic,
				'isWall' : isWall,
				'density' : density,
				'mapInvert': mapInvert,
				'mapThickness': mapThickness,
				'mapResolution': mapResolution,
				'collisionObjectType': 5,
				'collisionObjectScale': [1,1,1],
				'resolutionSDF': mapResolution,
				'invertSDF': mapInvert
			}
	 
		scene['RigidBodies'].append(rb)
		
	######################################################
	# add fluid
	######################################################
	def addFluid(self, scene, fluidNode, scenePath):
		
		# export geometry
		cmds.select([fluidNode], replace=True)
		polyTri = cmds.polyTriangulate()
		name = cmds.ls( selection=True, type='transform', long=False )[0] 
		fileName = os.path.join(scenePath, "fluid_" + name + ".obj")
		cmds.file(fileName, force=True, options="groups=0;ptgroups=0;materials=0;smoothing=0;normals=0", pr=True, exportSelected=True, type="OBJexport")
		cmds.delete(polyTri)

		initialVelocity = cmds.getAttr(fluidNode + ".SPH_initialVelocity")[0]
		resSDF = cmds.getAttr(fluidNode + ".SPH_resolutionSDF")[0]
		fluidId = cmds.getAttr(fluidNode + ".SPH_fluidId")
		
		# get all fluidIds
		fluidIds = cmds.ls( type='SPHFluidConfigurationNode', long=False ) 
		if fluidId not in fluidIds:
			cmds.warning("The node '" + fluidNode + "' uses the fluid material '" + fluidId + "' which is not defined in the scene.")

		denseMode = cmds.getAttr(fluidNode + ".SPH_denseMode")
		invert = cmds.getAttr(fluidNode + ".SPH_invert")
		fluid = {  
				'id': fluidId,
				'particleFile': "fluid_" + name + ".obj",
				'translation': [0,0,0],
				'rotationAxis': [1,0,0],
				'rotationAngle': 0.0,
				'scale': [1,1,1],
				'initialVelocity': initialVelocity,
				'resolutionSDF': resSDF,
				'denseMode': denseMode,
				'invert': invert
			}
	 
		scene['FluidModels'].append(fluid)
		
		

	######################################################
	# add rectangular emitter 
	######################################################
	def addRectangularEmitter(self, sphConfig, scene, node, scenePath):
		
		t = cmds.xform(node, query=True, t=True, ws=True)
		s = cmds.xform(node, query=True, s=True, ws=True)
		
		# get particleRadius
		radius = self.getCurrentParticleRadius(sphConfig)
		diam = 2.0 * radius
		s[1] -= 2.0*diam
		s[2] -= 2.0*diam
			
		axisAngle = PluginFunctions.getAxisAngle(node)
		
		startTime = cmds.getAttr(node + ".SPH_startTime")
		endTime = cmds.getAttr(node + ".SPH_endTime")
		velocity = cmds.getAttr(node + ".SPH_velocity")
		id = cmds.getAttr(node + ".SPH_fluidId")
		emitter = {  
				'id': id,
				'width': int(s[2]/diam), 
				'height': int(s[1]/diam),
				'translation': t,
				'rotationAxis': axisAngle[0],
				'rotationAngle': axisAngle[1],
				'startTime': startTime,
				'endTime': endTime,
				'velocity' : velocity,
				'type' : 0
			}
	 
		scene['Emitters'].append(emitter)
		

	######################################################
	# add circular emitter 
	######################################################
	def addCircularEmitter(self, sphConfig, scene, node, scenePath):
		
		t = cmds.xform(node, query=True, t=True, ws=True)
		s = cmds.xform(node, query=True, s=True, ws=True)
		
		# get particleRadius
		radius = self.getCurrentParticleRadius(sphConfig)
		s[1] -= 2.0*radius
		
		axisAngle = PluginFunctions.getAxisAngle(node)
		
		startTime = cmds.getAttr(node + ".SPH_startTime")
		endTime = cmds.getAttr(node + ".SPH_endTime")
		velocity = cmds.getAttr(node + ".SPH_velocity")
		id = cmds.getAttr(node + ".SPH_fluidId")
		emitter = {  
				'id': id,
				'width': int(s[1]/radius), 
				'translation': t,
				'rotationAxis': axisAngle[0],
				'rotationAngle': axisAngle[1],
				'startTime': startTime,
				'endTime': endTime,
				'velocity' : velocity,
				'type' : 1
			}
	 
		scene['Emitters'].append(emitter)	

	######################################################
	# add animation field 
	######################################################
	def addAnimationField(self, sphConfig, scene, node, scenePath):
		
		t = cmds.xform(node, query=True, t=True, ws=True)
		s = cmds.xform(node, query=True, s=True, ws=True)
		axisAngle = PluginFunctions.getAxisAngle(node)
		
		particleField = cmds.getAttr(node + ".SPH_particleField")
		shapeType = cmds.getAttr(node + ".SPH_shapeType")
		expression_x = cmds.getAttr(node + ".SPH_expressionX")
		expression_y = cmds.getAttr(node + ".SPH_expressionY")
		expression_z = cmds.getAttr(node + ".SPH_expressionZ")
		
		animField = {  
				'particleField': particleField,
				'translation': t,
				'rotationAxis': axisAngle[0],
				'rotationAngle': axisAngle[1],
				'scale': s,
				'shapeType': shapeType,
				'expression_x' : expression_x,
				'expression_y' : expression_y,
				'expression_z' : expression_z
			}
	 
		scene['AnimationFields'].append(animField)		
							

# Node definition
class SPHConfigurationNode(OpenMayaMPx.MPxLocatorNode):
	kPluginNodeId = OpenMaya.MTypeId(0x90000)
	kPluginNodeTypeName = "SPHConfigurationNode"
	 
	# class variables
	input = OpenMaya.MObject()
	dataAttr = OpenMaya.MObject()
	sphParameters = OrderedDict()
	
	def __init__(self):
		OpenMayaMPx.MPxLocatorNode.__init__(self)
			
	
	def postConstructor(self):
		OpenMayaMPx.MPxLocatorNode.postConstructor(self) 
		
	# initializer
	@staticmethod
	def initialize():
	
		SPHConfigurationNode.initParameters()
		
		# add attributes	
		for key in SPHConfigurationNode.sphParameters:
			params = SPHConfigurationNode.sphParameters[key]
			for param in params:
				paramType = param["type"]
				paramName = param["name"]
				paramLabel = param["label"]
				if paramType == "bool":
					attr = PluginFunctions.createBoolAttr("SPH_" + paramName, paramName, param["value"])
					SPHConfigurationNode.addAttribute( attr )						
				elif paramType == "float":
					attr = PluginFunctions.createFloatAttr("SPH_" + paramName, paramName, param["value"], param["min"], param["max"], param["fieldMin"], param["fieldMax"])
					SPHConfigurationNode.addAttribute( attr )					
				elif paramType == "int":
					attr = PluginFunctions.createIntAttr("SPH_" + paramName, paramName, param["value"], param["min"], param["max"], param["fieldMin"], param["fieldMax"])
					SPHConfigurationNode.addAttribute( attr )
				elif paramType == "vec3":
					attr = PluginFunctions.createVec3Attr("SPH_" + paramName, paramName, param["value"])
					SPHConfigurationNode.addAttribute( attr )
				elif paramType == "enum":
					attr = PluginFunctions.createEnumAttr("SPH_" + paramName, paramName, param["value"], param["enumList"])
					SPHConfigurationNode.addAttribute( attr )
		

	# creator
	@staticmethod
	def creator():
		return OpenMayaMPx.asMPxPtr( SPHConfigurationNode() )

	def compute(self,plug,dataBlock):
		# if ( plug == SPHConfigurationNode.output ):
				# dataHandle = dataBlock.inputValue( SPHConfigurationNode.input )
				
				# inputFloat = dataHandle.asFloat()
				# result = math.sin( inputFloat ) * 10.0
				# outputHandle = dataBlock.outputValue( SPHConfigurationNode.output )
				# outputHandle.setFloat( result )
				# dataBlock.setClean( plug )

		return OpenMaya.kUnknownParameter

	######################################################
	# initParameters
	######################################################	
	@staticmethod	
	def initParameters():
		SPHConfigurationNode.sphParameters["General"] = [
			PluginFunctions.createBoolParam("pause", "Pause", "Pause simulation after loading.", True),
			PluginFunctions.createFloatParam("timeStepSize", "Time step size", "Time step size", 0.001, 0.00001, 1.0),
			PluginFunctions.createFloatParam("pauseAt", "Pause simulation at", "Pause simulation at the given time. When the value is negative, the simulation is not paused.", -1, -1, 100, -1),
			PluginFunctions.createFloatParam("stopAt", "Stop simulation at", "Stop simulation at the given time. When the value is negative, the simulation is not stopped.", -1, -1, 100, -1)
		]
		SPHConfigurationNode.sphParameters["Visualization"] = [	
			PluginFunctions.createVec3Param("cameraPosition", "Camera position", "Initial position of the camera.", [0.0,3.0,8.0]),
			PluginFunctions.createVec3Param("cameraLookat", "Camera lookat", "Lookat point of the camera.", [0.0,0.0,0.0]),
			PluginFunctions.createIntParam("numberOfStepsPerRenderUpdate", "# time steps / update", "Number of simulation steps per rendered frame.", 4, 1, 100),
			PluginFunctions.createEnumParam("renderWalls", "Render walls", "Make walls visible/invisible.", 4, ["None", "Particles (all)", "Particles (no walls)", "Geometry (all)", "Geometry (no walls)"]),
		]
		SPHConfigurationNode.sphParameters["Export"] = [	
			PluginFunctions.createBoolParam("enablePartioExport", "Partio export", "Enable/disable partio export.", False),
			PluginFunctions.createBoolParam("enableVTKExport", "VTK export", "Enable/disable VTK export.", False),
			PluginFunctions.createFloatParam("dataExportFPS", "Export FPS", "Frame rate of particle export.", 25, 0.1, 1000),
			PluginFunctions.createStringParam("particleAttributes", "Export attributes", "Attributes that are exported in the particle files (except id and position).", "velocity"),
			PluginFunctions.createBoolParam("enableStateExport", "State export", "Enable/disable simulation state export.", False),
			PluginFunctions.createFloatParam("stateExportFPS", "State export FPS", "Frame rate of state export.", 1, 0.1, 1000)
		]
		SPHConfigurationNode.sphParameters["Simulation"] = [
			PluginFunctions.createBoolParam("sim2D", "2D simulation", "2D/3D simulation.", False),
			PluginFunctions.createBoolParam("enableZSort", "Enable z-sort", "Enable z-sort to improve cache hits.", True),
			PluginFunctions.createFloatParam("particleRadius", "Particle radius", "Radius of the fluid particles.", 0.025, 0.0001, 1000.0, 0),
			PluginFunctions.createVec3Param("gravitation", "Gravitation", "Vector to define the gravitational acceleration.", [0,-9.81,0]),
			PluginFunctions.createEnumParam("simulationMethod", "Simulation method", "Simulation method.", 4, ["WCSPH", "PCISPH", "PBF", "IISPH", "DFSPH", "Projective Fluids"]),
			PluginFunctions.createIntParam("maxIterations", "Max. iterations", "Maximal number of iterations of the pressure solver.", 100, 1, 1000, 1),
			PluginFunctions.createFloatParam("maxError", "Max. density error(%)", "Maximal density error (%).", 0.01, 1.0e-6, 1.0, 0),
			PluginFunctions.createEnumParam("boundaryHandlingMethod", "Boundary handling method", "Boundary handling method.", 2, ["Akinci et al. 2012", "Koschier and Bender 2017", "Bender et al. 2019"])
		
		]
		SPHConfigurationNode.sphParameters["CFL"] = [	
			PluginFunctions.createEnumParam("cflMethod", "CFL - method", "CFL method used for adaptive time stepping.", 1, ["None", "CFL", "CFL - iterations"]),
			PluginFunctions.createFloatParam("cflFactor", "CFL - factor", "Factor to scale the CFL time step size.", 0.5, 1e-6, 10.0, 0),
			PluginFunctions.createFloatParam("cflMaxTimeStepSize", "CFL - max. time step size", "Max. time step size.", 0.005, 1e-6, 1.0, 0)
		]
		SPHConfigurationNode.sphParameters["Kernel"] = [	
			PluginFunctions.createEnumParam("kernel", "Kernel", "Kernel function used in the SPH model (in 2D use only cubic or Wendland).", 4, ["Cubic spline", "Wendland quintic C2", "Poly6", "Spiky", "Precomputed cubic spline"]),
			PluginFunctions.createEnumParam("gradKernel", "Gradient of kernel", "Gradient of the kernel function used in the SPH model (in 2D use only cubic or Wendland).", 4, ["Cubic spline", "Wendland quintic C2", "Poly6", "Spiky", "Precomputed cubic spline"])
		]
		SPHConfigurationNode.sphParameters["WCSPH"] = [
			PluginFunctions.createFloatParam("stiffness", "Stiffness", "Stiffness coefficient of EOS.", 10000, 0, 500000),
			PluginFunctions.createFloatParam("exponent", "Exponent (gamma)", "Exponent of EOS.", 7.0, 1.0e-6, 10.0, 0)		
		]
		SPHConfigurationNode.sphParameters["PBF"] = [
			PluginFunctions.createEnumParam("velocityUpdateMethod", "Velocity update method", "Method for the velocity integration.", 0, ["First Order Update", "Second Order Update"])	
		]
		SPHConfigurationNode.sphParameters["DFSPH"] = [
			PluginFunctions.createIntParam("maxIterationsV", "Max. iterations (divergence)", "Maximal number of iterations of the divergence solver.", 100, 1, 1000, 1),
			PluginFunctions.createFloatParam("maxErrorV", "Max. divergence error(%)", "Maximal divergence error (%).", 0.01, 1.0e-6, 1.0, 0),
			PluginFunctions.createBoolParam("enableDivergenceSolver", "Enable divergence solver", "Turn divergence solver on/off.", True)		
		]
		SPHConfigurationNode.sphParameters["Projective Fluids"] = [
			PluginFunctions.createFloatParam("stiffnessPF", "Stiffness", "Stiffness coefficient.", 50000, 0, 500000)
		]

# Node definition
class SPHFluidConfigurationNode(OpenMayaMPx.MPxLocatorNode):
	kPluginNodeId = OpenMaya.MTypeId(0x90001)
	kPluginNodeTypeName = "SPHFluidConfigurationNode"
	 
	# class variables
	input = OpenMaya.MObject()
	dataAttr = OpenMaya.MObject()
	sphParameters = OrderedDict()
	
	def __init__(self):
		OpenMayaMPx.MPxLocatorNode.__init__(self)
				
	# initializer
	@staticmethod
	def initialize():
	
		SPHFluidConfigurationNode.initParameters()
		
		# add attributes	
		for key in SPHFluidConfigurationNode.sphParameters:
			params = SPHFluidConfigurationNode.sphParameters[key]
			for param in params:
				paramType = param["type"]
				paramName = param["name"]
				paramLabel = param["label"]
				if paramType == "bool":
					attr = PluginFunctions.createBoolAttr("SPH_" + paramName, paramName, param["value"])
					SPHConfigurationNode.addAttribute( attr )						
				elif paramType == "float":
					attr = PluginFunctions.createFloatAttr("SPH_" + paramName, paramName, param["value"], param["min"], param["max"], param["fieldMin"], param["fieldMax"])
					SPHConfigurationNode.addAttribute( attr )					
				elif paramType == "int":
					attr = PluginFunctions.createIntAttr("SPH_" + paramName, paramName, param["value"], param["min"], param["max"], param["fieldMin"], param["fieldMax"])
					SPHConfigurationNode.addAttribute( attr )
				elif paramType == "vec3":
					attr = PluginFunctions.createVec3Attr("SPH_" + paramName, paramName, param["value"])
					SPHConfigurationNode.addAttribute( attr )
				elif paramType == "enum":
					attr = PluginFunctions.createEnumAttr("SPH_" + paramName, paramName, param["value"], param["enumList"])
					SPHConfigurationNode.addAttribute( attr )
		

	# creator
	@staticmethod
	def creator():
		return OpenMayaMPx.asMPxPtr( SPHFluidConfigurationNode() )

	######################################################
	# initParameters
	######################################################	
	@staticmethod	
	def initParameters():
		SPHFluidConfigurationNode.sphParameters["Simulation"] = [
			PluginFunctions.createFloatParam("density0", "Rest density", "Rest density of the fluid.", 1000.0, 0.1, 10000.0)
		]	
		SPHFluidConfigurationNode.sphParameters["Visualization"] = [	
			PluginFunctions.createStringParam("colorField", "Color field", "Choose vector or scalar field for particle coloring.", "velocity"),
			PluginFunctions.createEnumParam("colorMapType", "Color map type", "Selection of a color map for coloring the scalar/vector field.", 1, ["None", "Jet", "Plasma"]),
			PluginFunctions.createFloatParam("renderMinValue", "Min. value", "Minimal value used for color-coding the color field in the rendering process.", 0, -1000, 1000, -1000000),
			PluginFunctions.createFloatParam("renderMaxValue", "Max. value", "Maximal value used for color-coding the color field in the rendering process.", 5, -1000, 1000, -1000000)
		]	
		SPHFluidConfigurationNode.sphParameters["Emitters"] = [	
			PluginFunctions.createIntParam("maxEmitterParticles", "Max. number of emitted particles", "Maximum number of emitted particles", 10000, 1, 10000000),
			PluginFunctions.createBoolParam("emitterReuseParticles", "Reuse particles", "Reuse particles if they are outside of the bounding box defined by emitterBoxMin, emitterBoxMaRex.", False),
			PluginFunctions.createVec3Param("emitterBoxMin", "Emitter box min.", "Minimum coordinates of an axis-aligned box (used in combination with emitterReuseParticles).", [0.0,0.0,0.0]),
			PluginFunctions.createVec3Param("emitterBoxMax", "Emitter box max.", "Maximum coordinates of an axis-aligned box (used in combination with emitterReuseParticles).", [1.0,1.0,1.0])
		]
		SPHFluidConfigurationNode.sphParameters["Viscosity"] = [	
			PluginFunctions.createEnumParam("viscosityMethod", "Viscosity", "Method to compute viscosity forces.", 1, ["None", "Standard", "XSPH", "Bender and Koschier 2017", "Peer et al. 2015", "Peer et al. 2016", "Takahashi et al. 2015 (improved)", "Weiler et al. 2018"]),
			PluginFunctions.createFloatParam("viscosity", "Viscosity coefficient", "Coefficient for the viscosity force computation.", 0.01, 0, 1000, 0),
			PluginFunctions.createIntParam("viscoMaxIter", "Max. iterations (visco)", "(Implicit solvers) Max. iterations of the viscosity solver.", 100, 1, 1000),
			PluginFunctions.createFloatParam("viscoMaxError", "Max. visco error", "(Implicit solvers) Max. error of the viscosity solver.", 0.01, 1e-6, 1, 0),
			PluginFunctions.createIntParam("viscoMaxIterOmega", "Max. iterations (vorticity diffusion)", "(Peer et al. 2016) Max. iterations of the vorticity diffusion solver.", 100, 1, 1000),
			PluginFunctions.createFloatParam("viscoMaxErrorOmega", "Max. vorticity diffusion error", "(Peer et al. 2016) Max. error of the vorticity diffusion solver.", 0.01, 1e-6, 1, 0),
			PluginFunctions.createFloatParam("viscosityBoundary", "Viscosity coefficient (Boundary)", "Coefficient for the viscosity force computation at the boundary.", 0.0, 0, 1000, 0)
		]
		SPHFluidConfigurationNode.sphParameters["Vorticity"] = [	
			PluginFunctions.createEnumParam("vorticityMethod", "Vorticity method", "Method to compute vorticity forces.", 0, ["None", "Micropolar model", "Vorticity confinement"]),
			PluginFunctions.createFloatParam("vorticity", "Vorticity coefficient", "Coefficient for the vorticity force computation.", 0.01, 0, 10.0, 0),
			PluginFunctions.createFloatParam("viscosityOmega", "Angular viscosity coefficient", "Viscosity coefficient for the angular velocity field.", 0.1, 0, 10.0, 0),
			PluginFunctions.createFloatParam("inertiaInverse", "Inertia inverse", "Inverse microinertia used in the micropolar model.", 0.5, 0, 10.0, 0)
		]
		SPHFluidConfigurationNode.sphParameters["Drag force"] = [	
			PluginFunctions.createEnumParam("dragMethod", "Drag method", "Method to compute drag forces.", 0, ["None", "Macklin et al. 2014", "Gissler et al. 2017"]),
			PluginFunctions.createFloatParam("drag", "Drag coefficient", "Coefficient for the drag force computation.", 0.01, 0, 100.0, 0)
		]
		SPHFluidConfigurationNode.sphParameters["Surface tension"] = [	
			PluginFunctions.createEnumParam("surfaceTensionMethod", "Surface tension method", "Method to compute surface tension forces.", 0, ["None", "Becker & Teschner 2007", "Akinci et al. 2013", "He et al. 2014"]),
			PluginFunctions.createFloatParam("surfaceTension", "Surface tension coefficient", "Coefficient for the surface tension computation.", 0.05, 0, 100.0, 0)
		]

######################################################
# loadRigidBodies
#
# load rigid body data that was exported by 
# a SPH simulation
######################################################		
class loadRigidBodiesCmd(OpenMayaMPx.MPxCommand):
	s_name = "loadRigidBodies"

	def __init__(self):
		OpenMayaMPx.MPxCommand.__init__(self)

	@staticmethod
	def creator():
		return loadRigidBodiesCmd()

	def doIt(self, args):
		self.addedNodes = []
		self.firstFileName = cmds.fileDialog2(ff="*.bin", fm=1, dir="")[0]
		indexlist = re.findall(r'\d+', self.firstFileName)
		if len(indexlist) == 0:
			cmds.warning("No frame index found in file name.")
			return
		
		self.firstFrame = int(indexlist[-1])
		self.redoIt()


	def redoIt(self):	
		self.loadRigidBodies()
	
	def loadRigidBodies(self):
	
		folderName = os.path.dirname(self.firstFileName)
		frameNumber = self.firstFrame
		firstFile = open(self.firstFileName, 'rb')
		
		# read number of bodies
		bytes = firstFile.read()
		firstFile.close()
		
		(numBodies,), bytes = struct.unpack('i', bytes[:4]), bytes[4:]
		
		objFiles = []
		transformNodes = []
		for i in range(0, numBodies):
			# determine length of file name string
			(strLength,), bytes = struct.unpack('i', bytes[:4]), bytes[4:]

			# read file name
			objFile, bytes = bytes[:strLength], bytes[strLength:]
			
			# Check for duplicates and create instances
			if objFile in objFiles:
				idx = objFiles.index(objFile)
				newNodes = cmds.duplicate(transformNodes[idx], instanceLeaf= True)
				transformNodes.append(newNodes[0])
				self.addedNodes.append(newNodes)
			else:
				objFileName = os.path.join(folderName, objFile)
				newNodes = cmds.file(objFileName, i=True, rnn=True, type="OBJ", options="mo=1") 
				
				transformNodes.append(newNodes[0])
				objFiles.append(objFile)
				self.addedNodes.append(newNodes)
					
				
			# Read scaling factors in first file
			(sx,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
			(sy,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
			(sz,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
			
			cmds.scale(sx, sy, sz, transformNodes[i])
				
			(isWall,), bytes = struct.unpack('?', bytes[:1]), bytes[1:]
			(colr,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
			(colg,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
			(colb,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
			(cola,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
			
			if isWall:
				cmds.setAttr((transformNodes[i] + '.visibility'), 0)
			
			cmds.setKeyframe(transformNodes[i], at="s", t=1)
			if frameNumber > 1:
				cmds.setKeyframe(transformNodes[i], at="visibility", t=1, value=0)	
				if not isWall:
					cmds.setKeyframe(transformNodes[i], at="visibility", t=frameNumber, value=1)
				
		
		# load transformations
		for i in range(0, numBodies):
			# Read translation in first file
			(x,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
			(y,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
			(z,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]

			# Read rotation in first file
			r = []
			for j in range(0,9):
				(value,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
				r.append(value)

			cmds.xform(transformNodes[i], p=True, m=[r[0],r[1],r[2],0,r[3],r[4],r[5],0,r[6],r[7],r[8],0,x,y,z,1])
			cmds.setKeyframe(transformNodes[i], at="t", t=frameNumber)
			cmds.setKeyframe(transformNodes[i], at="r", t=frameNumber)
		
		# read other files
		idx = self.firstFileName.rfind(str(frameNumber))
		l = len(str(frameNumber))
		chk = True
		while chk:
			frameNumber += 1
			fileName = str(self.firstFileName[0:idx]) + str(frameNumber) + str(self.firstFileName[idx+l:])
			chk = os.path.exists(fileName)
			if chk:
				f = open(fileName, 'rb')
				bytes = f.read()
				f.close()
				
				# load transformations
				for i in range(0, numBodies):
					# Read translation in file
					(x,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
					(y,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
					(z,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]

					# Read rotation in file
					r = []
					for j in range(0,9):
						(value,), bytes = struct.unpack('f', bytes[:4]), bytes[4:]
						r.append(value)

					cmds.xform(transformNodes[i], p=True, m=[r[0],r[1],r[2],0,r[3],r[4],r[5],0,r[6],r[7],r[8],0,x,y,z,1])
					cmds.setKeyframe(transformNodes[i], at="t", t=frameNumber)
					cmds.setKeyframe(transformNodes[i], at="r", t=frameNumber)
					
		cmds.currentTime(1)
		
	def undoIt(self):
		for node in self.addedNodes:
			print node
			cmds.delete(node)

	def isUndoable(self):
		return True

         
######################################################
# createSPHMenu
######################################################		
def createSPHMenu():
	global menuId
	menuId = cmds.menu( label='SPlisHSPlasH', p="MayaWindow" )
	cmds.menuItem(divider=True, dividerLabel="Scene generating")

	cmds.menuItem( label='Add scene configuration',command=
		'if "SPH_Config" not in cmds.ls( type="transform"):\n' +
		'	cmds.createNode("transform", name="SPH_Config")\n' +
		'cmds.createNode("SPHConfigurationNode", name="Configuration", parent="SPH_Config")')
		
	cmds.menuItem( label='Add fluid material',command=
		'if "SPH_Fluid" not in cmds.ls( type="transform"):\n' +
		'	cmds.createNode("transform", name="SPH_Fluid")\n' +
		'cmds.createNode("SPHFluidConfigurationNode", name="Fluid", parent="SPH_Fluid")')
		
	cmds.menuItem(divider=True)
	cmds.menuItem( label='Convert selection to fluid',command='cmds.convertToFluid()' )
	cmds.menuItem( label='Convert selection to rigid bodies',command='cmds.convertToRigidBodies()' )
	cmds.menuItem(divider=True)
	cmds.menuItem( label='Create rectangular emitter',command='cmds.createRectangularEmitter()' )
	cmds.menuItem( label='Create circular emitter',command='cmds.createCircularEmitter()' )
	cmds.menuItem(divider=True)
	cmds.menuItem( label='Create box animation field',command='cmds.createAnimationField(s=0)' )
	cmds.menuItem( label='Create sphere animation field',command='cmds.createAnimationField(s=1)' )
	cmds.menuItem( label='Create cylinder animation field',command='cmds.createAnimationField(s=2)' )
	cmds.menuItem(divider=True)
	cmds.menuItem( label='Save scene',command='cmds.saveModel()' )
	cmds.menuItem(divider=True, dividerLabel="Import")
	cmds.menuItem( label='Load rigid body data',command='cmds.loadRigidBodies()' )


######################################################
# deleteSPHMenu
######################################################	
def deleteSPHMenu():
	global menuId
	cmds.deleteUI(menuId)
	return
	

# Initialize the script plug-in
def initializePlugin(mobject):
	global settingsWinId
	global fluidWinId
	global menuId
	global fluidIds
	global sphParameters
	global fluidParameters
	
	mplugin = OpenMayaMPx.MFnPlugin(mobject, "SPlisHSPlasH", "1.0", "Any")

	settingsWinId = ""
	fluidWinId = ""
	menuId = ""
	fluidIds = ["Fluid"]
 	
	try:
		mplugin.registerNode( SPHConfigurationNode.kPluginNodeTypeName, SPHConfigurationNode.kPluginNodeId, SPHConfigurationNode.creator, SPHConfigurationNode.initialize, OpenMayaMPx.MPxNode.kLocatorNode )
		mplugin.registerNode( SPHFluidConfigurationNode.kPluginNodeTypeName, SPHFluidConfigurationNode.kPluginNodeId, SPHFluidConfigurationNode.creator, SPHFluidConfigurationNode.initialize, OpenMayaMPx.MPxNode.kLocatorNode )
		
		mplugin.registerCommand(createRectangularEmitterCmd.s_name, createRectangularEmitterCmd.creator)
		mplugin.registerCommand(createCircularEmitterCmd.s_name, createCircularEmitterCmd.creator)
		mplugin.registerCommand(saveModelCmd.s_name, saveModelCmd.creator)
		mplugin.registerCommand(convertToFluidCmd.s_name, convertToFluidCmd.creator)
		mplugin.registerCommand(convertToRigidBodiesCmd.s_name, convertToRigidBodiesCmd.creator)
		mplugin.registerCommand(createAnimationFieldCmd.s_name, createAnimationFieldCmd.creator, createAnimationFieldCmd.syntaxCreator)	
		mplugin.registerCommand(loadRigidBodiesCmd.s_name, loadRigidBodiesCmd.creator)	
	
		
	except:
		sys.stderr.write( "Failed to register nodes." )
		raise
		
	createSPHMenu()

# Uninitialize the script plug-in
def uninitializePlugin(mobject):
	mplugin = OpenMayaMPx.MFnPlugin(mobject)
	deleteSPHMenu()
	
	try:	
		mplugin.deregisterCommand(createRectangularEmitterCmd.s_name)
		mplugin.deregisterCommand(createCircularEmitterCmd.s_name)
		mplugin.deregisterCommand(saveModelCmd.s_name)
		mplugin.deregisterCommand(convertToFluidCmd.s_name)
		mplugin.deregisterCommand(convertToRigidBodiesCmd.s_name)
		mplugin.deregisterCommand(createAnimationFieldCmd.s_name)
		mplugin.deregisterCommand(loadRigidBodiesCmd.s_name)
		
		mplugin.deregisterNode( SPHFluidConfigurationNode.kPluginNodeId )
		mplugin.deregisterNode( SPHConfigurationNode.kPluginNodeId )
	except:
		sys.stderr.write( "Failed to deregister node")
		raise


