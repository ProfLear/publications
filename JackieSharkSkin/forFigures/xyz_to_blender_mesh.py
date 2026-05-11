#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 06:00:19 2023

@author: benjaminlear
"""
#% Blender stuff...
import bpy
from pathlib import Path
import math

folder = Path(r"/Users/benjaminlear/GitHub/publications/JackieSharkSkin/forFigures")
#stem = Path("smallPreBlender")

def make_scale_bar( x, y, z, length, height, name = "scale bar"):
    if f"{name} scale bar" not in bpy.data.objects:
        bpy.ops.mesh.primitive_plane_add(size=1, enter_editmode=False, align='WORLD', location=(x, y, z))
        plane = bpy.context.active_object
        #
        # Link the created object to the surface collection
        surface_collection.objects.link(plane)
        bpy.context.scene.collection.objects.unlink(plane)
        #
        # Rename the object to "scale bar"
        plane.name = f"{name} scale bar"
        #
        # Scale the plane in x, y, z dimensions
        bpy.data.objects[f"{name} scale bar"].scale = (length, height, 1)
        #
        # rotate the plane
        bpy.data.objects[f"{name} scale bar"].rotation_euler[0] = math.radians(90)  # Rotation around X-axis
        #
        # Update the plane's transformation
        bpy.data.objects[f"{name} scale bar"].update_from_editmode()
        #
        # Create a new material with a pure red color
        material_red = bpy.data.materials.new(name="Scale Bar")
        material_red.diffuse_color = (1, 0, 0, 1)  # RGBA, where (1, 0, 0) is red
        #
        # Assign the material to the plane
        if bpy.data.objects[f"{name} scale bar"].data.materials:
            # Assign to 1st material slot
            bpy.data.objects[f"{name} scale bar"].data.materials[0] = material_red
        else:
            # No slots
            bpy.data.objects[f"{name} scale bar"].data.materials.append(material_red)
        #
        # make sure it is visible in render
        bpy.data.objects[f"{name} scale bar"].hide_render = False
        #
        # TEXT LABEL
        #
        # Create a new text object
        bpy.ops.object.text_add(location=(0, 0, 0))
        label = bpy.context.object
        label.data.body = "10 microns"
        label.name = f"{name} label"
        # Link the created object to the surface collection
        surface_collection.objects.link(label)
        bpy.context.scene.collection.objects.unlink(label)
        #
        # Set alignment to 'CENTER' for horizontal and 'TOP' for vertical
        label.data.align_x = 'CENTER'
        label.data.align_y = 'TOP_BASELINE'
        #
        # Adjust the location, rotation, and scale as needed
        label.location = (x, y, z - 0.2) # Set location
        label.rotation_euler = (math.radians(90), 0, 0)  # Set rotation (in radians)
        label.scale = (0.2, 0.2, 1)  # Set scale
        #
        if label.data.materials:
            # Assign to 1st material slot
            label.data.materials[0] = material_red
        else:
            # No slots
            label.data.materials.append(material_red)
    


def verts_faces_from_csv(file, meshname = "Mesh"):
    verts = []
    faces= []
    x_min = 0
    x_max = 0
    y_min = 0
    y_max = 0
    z_min = 0
    with open(file, "r") as f:
        for row in f:
            readrow = row.split(",")
            if len(readrow) == 4:
                temprow = []
                for e in readrow[0:3]:
                    temprow.append(float(e)/100) # scale so that each meter is 100 microns
                if temprow[0] < x_min:
                    x_min = temprow[0]
                if temprow[0] > x_max:
                    x_max = temprow[0]
                if temprow[1] < y_min:
                    y_min = temprow[1]
                if temprow[1] > y_min:
                    y_max = temprow[1]
                if temprow[2] < z_min:
                    z_min = temprow[2]
                temprow[2]  = temprow[2]*100 # change this back for the z-axis so that 1 meter is 1 micron
                verts.append(temprow)
            if len(readrow) == 5:
                temprow = []
                for e in readrow[0:4]:
                    temprow.append(int(e))
                faces.append(temprow)
        #        
        # now let us go and normalize the x-y data to the x-axis, this will make the x-axis 10 meters long,and center it at x=0
        x_range = x_max - x_min
        y_range = y_max - y_min
        for i, row in enumerate(verts):
            verts[i][0] = ((verts[i][0] - x_min)/x_range - 0.5) * 10 # normalize, center on x=0 and scale
            verts[i][1] = (verts[i][1]/x_range - 0.5 * y_range/x_range) * 10 # normalize wrt x_range, center on y=0 and scale
            verts[i][2] = (verts[i][2])/x_range * 10 # already centered as we want it... so, just normalize wrt the x-range
    #Delete the default objects
    #add the sun
    # add scale bar
    #print(z_min)
    make_scale_bar(-5 + 1/x_range*10/2, # x position
                   -5 - 0.05, # y position
                   z_min/x_range*10*100 - 0.05, # z position
                   1/x_range * 10, # length
                   0.1/x_range*10, # height
                   name = meshname)
    #      
    # Create a mesh and an object
    mesh = bpy.data.meshes.new(meshname)
    obj = bpy.data.objects.new(meshname, mesh)
    #
    # Link the created mesh object to the surface collection
    mesh_obj = bpy.data.objects[meshname]
    surface_collection.objects.link(mesh_obj)
    #
    # Create the mesh from the vertices and faces
    mesh.from_pydata(verts, [], faces)
    #
    # Update the mesh with the new data
    mesh.update()     


def add_sun(name="Sun", power = 1, rotation=(0, 0, 0)):
    # Create a new sun lamp
    bpy.ops.object.light_add(type='SUN')
    sun = bpy.context.object
    #
    # Set the name of the sun
    sun.name = name
    #
    # Set the sun's energy and location
    sun.data.energy = power  # Adjust the energy level as needed
    sun.location = (0, 0, 1)  # Adjust the location as needed
    #
    # Set the sun's rotation (Euler angles in radians)
    sun.rotation_euler[0] = math.radians(rotation[0])  # Rotation around X-axis
    sun.rotation_euler[1] = math.radians(rotation[1])  # Rotation around Y-axis
    sun.rotation_euler[2] = math.radians(rotation[2])  # Rotation around Z-axis
    
def add_camera(name = "Camera", location = (0, -19, 5), rotation = (75, 0 , 0)):
    # Create a new camera
    bpy.ops.object.camera_add()
    camera = bpy.context.object
    #
    camera.name = name
    # Set the camera's location and rotation
    camera.location = location  # Adjust the location as needed
    camera.rotation_euler[0] = math.radians(rotation[0])
    camera.rotation_euler[1] = math.radians(rotation[1])
    camera.rotation_euler[2] = math.radians(rotation[2])
    #
    # Link the camera to the scene
    scene = bpy.context.scene
    #bpy.context.scene.collection.objects.link(camera)
    #
    # Set the newly added camera as the active camera
    scene.camera = camera


# Delete default objects at the start
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete()

# Check if the default collection is the master collection
default_collection = bpy.data.collections.get("Collection")
if default_collection and default_collection != bpy.context.scene.collection:
    # Remove all objects from the default collection
    for obj in default_collection.objects:
        default_collection.objects.unlink(obj)
    #
    # Delete the default collection
    bpy.data.collections.remove(default_collection)
       

for surface in ["smallPre", "largePre", "smallPost", "largePost", "flatPre", "flatPost"]:
    # Create a new collection for the surface
    surface_collection = bpy.data.collections.new(surface)
    bpy.context.scene.collection.children.link(surface_collection)
    #
    # Read vertices and faces from CSV and create the mesh
    meshname = surface
    verts_faces_from_csv(folder/f"{surface}.csv", meshname=meshname)


add_sun(name = "main", power = 8, rotation = (38, -83, 115))
add_sun(name = "auxillary", power = 0.5, rotation = (0, 90, 225))
add_sun(name = "backlight1", power = 0.5, rotation = (0, 90, 45))
add_sun(name = "backlight2", power = 0.5, rotation = (0, 90, 135))

add_camera(name = "25 degrees")

# Access the current scene
scene = bpy.context.scene

# Set render settings for transparency and to cycles
scene.render.film_transparent = True
bpy.context.scene.render.engine = 'CYCLES'



