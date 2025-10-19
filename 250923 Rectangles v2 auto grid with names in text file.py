"""
# -*- coding: utf-8 -*-
Rhino Roof Generator - Creates architectural roof structure with walls and dual doors
Uses proven geometric methods with organized formatting
"""

import rhinoscriptsyntax as rs
import Rhino
import scriptcontext as sc
import math
import os
import subprocess

def log_specs_to_file(text, filepath="C:/temp/building_specs_log.txt", auto_open=False):
    """Append the specs text to a log file."""
    # Ensure folder exists
    folder = os.path.dirname(filepath)
    if not os.path.exists(folder):
        os.makedirs(folder)
    try:
        with open(filepath, 'a') as f:   # 'a' = append
            f.write(text + "\n\n")       # blank line between entries
    except Exception as e:
        print("Could not write log:", e)
    if auto_open:
        try:
            # This will open with system’s default text editor (Notepad)
            os.startfile(filepath)
        except:
            # fallback for older setups
            subprocess.Popen(['notepad.exe', filepath])

# =============================================================================
# CONFIGURATION - Proven Working Parameters
# =============================================================================

PARAMS = {
    # === BASE RECTANGLE SETTINGS ===
    "center": (0.0, 0.0, 0.0),    # Building center point (x, y, z)
    "length": 40.0,               # Building length along X-axis
    "width": 24.0,                # Building width along Y-axis
    
    # === RIDGE LINE CONFIGURATION ===
    # Controls where the roof peak line is positioned
    "edge_index": 0,              # Which edge to start ridge from (0-3)
    "t0": 0.7,                    # Position along first edge (0.0-1.0)
    "t1": 0.75,                   # Position along opposite edge (0.0-1.0)
    
    # === ROOF ANGLES & ELEVATION ===
    # Creates the characteristic gable roof shape
    "rotate_deg_pos": -25.0,      # Rotation angle for positive side (degrees)
    "rotate_deg_neg": 10.0,       # Rotation angle for negative side (degrees)
    "ridge_height": 20.0,         # Final roof elevation above ground
    
    # === WALL INSET DISTANCES ===
    # Creates interior walls smaller than the base rectangle
    "offset_A": 10.0,             # Left side inset distance
    "offset_B": 7.0,              # Right side inset distance
    "offset_C": 4.0,              # Bottom side inset distance
    "offset_D": 5.0,              # Top side inset distance
    
    # === WALL HEIGHT ===
    "b_extrude": 50.0,            # Vertical extrusion height for walls
    
    # === DUAL DOOR CONFIGURATION ===
    # Creates two doors that can be placed anywhere along building perimeter
    "ap_u": 0.25,                 # First door position (0.0-1.0)
    "ap_width": 6.0,              # First door width
    "ap_height": 7.0,             # First door height
    "ap_depth": 6.0,              # First door depth
    
    "ap2_u": 0.75,                # Second door position (opposite side)
    "ap2_width": 6.0,             # Second door width  
    "ap2_height": 7.0,            # Second door height
    "ap2_depth": 6.0,             # Second door depth
    
    # === SKYLIGHT CONFIGURATION ===
    # Creates skylight that cuts through roof (and walls if overlapping)
    "skylight_x": 5.0,            # X offset from building center
    "skylight_y": 5.0,            # Y offset from building center  
    "skylight_length": 8.0,       # Skylight length
    "skylight_width": 15.0,       # Skylight width
    "skylight_height": 3.0,      # Skylight cutting box height (penetrates roof and walls)
    "skylight_z": 0.0,            # Z offset from building center (for independent positioning)
    
    # === ADDITIONAL SURFACE (Rectangle C) ===
    "C_center": (15.0, -5.0, 0.0),
    "C_length": 20.0,
    "C_width": 15.0
}

# =============================================================================
# UTILITY FUNCTIONS - Core Mathematical Operations
# =============================================================================

def clamp01(x):
    """Clamp value to range [0.0, 1.0]"""
    return max(0.0, min(1.0, float(x)))

def point_on_curve_norm(crv_id, t01):
    """Get point on curve at normalized parameter"""
    crv = rs.coercecurve(crv_id)
    dom = crv.Domain
    u = dom.T0 + clamp01(t01) * (dom.T1 - dom.T0)
    return crv.PointAt(u)

def zmax_of(obj_id):
    """Get maximum Z coordinate of object's bounding box"""
    bb = rs.BoundingBox(obj_id)
    if not bb: 
        return float("inf")
    return max(p.Z for p in bb)

def bbox_center(bb):
    """Calculate center point of bounding box"""
    if not bb:
        return None
    x = sum(p.X for p in bb)/len(bb)
    y = sum(p.Y for p in bb)/len(bb)
    z = sum(p.Z for p in bb)/len(bb)
    return Rhino.Geometry.Point3d(x,y,z)

# =============================================================================
# GEOMETRY CREATION FUNCTIONS - Rhino Object Generation
# =============================================================================

def add_centered_rectangle(center, L, W, name=None):
    """Create rectangle centered at point"""
    cx, cy, cz = center
    plane = Rhino.Geometry.Plane(
        Rhino.Geometry.Point3d(cx, cy, cz),
        Rhino.Geometry.Vector3d.XAxis,
        Rhino.Geometry.Vector3d.YAxis
    )
    rx = Rhino.Geometry.Interval(-L*0.5, L*0.5)
    ry = Rhino.Geometry.Interval(-W*0.5, W*0.5)
    rect = Rhino.Geometry.Rectangle3d(plane, rx, ry)
    rid = sc.doc.Objects.AddCurve(rect.ToNurbsCurve())
    if name: 
        rs.ObjectName(rid, name)
    return rid

def add_rect_by_bounds(xmin, xmax, ymin, ymax, z=0.0, name=None):
    """Create rectangle by boundary coordinates"""
    pts = [(xmin,ymin,z), (xmax,ymin,z), (xmax,ymax,z), (xmin,ymax,z), (xmin,ymin,z)]
    gid = rs.AddPolyline(pts)
    if name: 
        rs.ObjectName(gid, name)
    return gid

def make_aperture_box_on_perimeter(boundary_id, u_norm, width, height, depth):
    """Create aperture box at any point along the joined rectangle perimeter.
    boundary_id: closed polyline curve id of rectangle boundary
    u_norm: normalized 0..1 parameter around full perimeter
    """
    crv = rs.coercecurve(boundary_id)
    if not crv:
        raise ValueError("Invalid boundary curve for aperture")

    total_len = crv.GetLength()
    s = clamp01(u_norm) * total_len

    ok, t = crv.LengthParameter(s)
    if not ok:
        dom = crv.Domain
        t = dom.T0 + (dom.T1-dom.T0)*clamp01(u_norm)

    pt = crv.PointAt(t)
    tan = crv.TangentAt(t)
    if tan.IsTiny():
        t2 = min(crv.Domain.T1, t + 1e-6*(crv.Domain.T1-crv.Domain.T0))
        tan = crv.TangentAt(t2)
    tan.Unitize()
    nrm = Rhino.Geometry.Vector3d(-tan.Y, tan.X, 0.0)

    half_w = width*0.5
    p0 = pt + tan*half_w
    p1 = pt - tan*half_w
    p2 = p1 + Rhino.Geometry.Vector3d(0,0,height)
    p3 = p0 + Rhino.Geometry.Vector3d(0,0,height)
    base_crv = rs.AddPolyline([p0, p1, p2, p3, p0])

    # Start at exterior wall surface, extrude inward only
    start_point = p0
    end_point = p0 + nrm*depth
    box_id = rs.ExtrudeCurveStraight(base_crv, start_point, end_point)
    
    if base_crv:
        rs.DeleteObject(base_crv)
    bb = rs.BoundingBox(box_id)
    return box_id, bb

def create_gable_roof(rectA_id, P):
    """Create gable roof using proven geometric method"""
    # Explode rectangle into edges
    edges = [e for e in rs.ExplodeCurves(rectA_id, True) or [] if e]
    
    # Get ridge points from opposite edges
    i0 = int(P["edge_index"]) % 4
    i1 = (i0 + 2) % 4
    pt0 = point_on_curve_norm(edges[i0], P["t0"])
    pt1 = point_on_curve_norm(edges[i1], P["t1"])
    
    # Create roof surface
    boundary = rs.JoinCurves(edges, True)[0]
    srf_id = rs.AddPlanarSrf(boundary)[0]
    
    # Create splitting plane (proven geometric method)
    v_axis = pt1 - pt0
    plane = Rhino.Geometry.Plane(pt0, v_axis, Rhino.Geometry.Vector3d(0,0,1))
    L, W = float(P["length"]), float(P["width"])
    R = max(L, W) * 2.0
    cutter_srf = Rhino.Geometry.PlaneSurface(
        plane, 
        Rhino.Geometry.Interval(-R,R), 
        Rhino.Geometry.Interval(-R,R)
    )
    cutter_id = sc.doc.Objects.AddSurface(cutter_srf)
    
    # Split roof surface
    halves = rs.SplitBrep(srf_id, cutter_id, True)
    rs.DeleteObject(cutter_id)
    
    if not halves or len(halves) != 2:
        print("Warning: Roof split failed")
        return srf_id
    
    # Rotate halves (tested geometric method)
    mid = Rhino.Geometry.Point3d(
        (pt0.X+pt1.X)/2.0, 
        (pt0.Y+pt1.Y)/2.0, 
        (pt0.Z+pt1.Z)/2.0
    )
    perp = Rhino.Geometry.Vector3d(-v_axis.Y, v_axis.X, 0.0)
    xform_pos = Rhino.Geometry.Transform.Rotation(
        math.radians(P["rotate_deg_pos"]), v_axis, pt0
    )
    xform_neg = Rhino.Geometry.Transform.Rotation(
        math.radians(P["rotate_deg_neg"]), v_axis, pt0
    )
    
    # Apply rotations based on which side each half is on
    for hid in halves:
        amp = Rhino.Geometry.AreaMassProperties.Compute(rs.coercebrep(hid))
        side = Rhino.Geometry.Vector3d.Multiply(amp.Centroid-mid, perp)
        sc.doc.Objects.Transform(hid, xform_pos if side>0 else xform_neg, True)
    
    # Join the rotated halves
    tol = sc.doc.ModelAbsoluteTolerance
    joined = Rhino.Geometry.Brep.JoinBreps(
        [rs.coercebrep(h) for h in halves], tol
    )
    roof_id = sc.doc.Objects.AddBrep(joined[0])
    rs.ObjectName(roof_id, "GableRoof")
    
    # Clean up
    for h in halves: 
        rs.DeleteObject(h)
    
    return roof_id

def create_walls_under_roof(P, roof_id):
    """Create walls that fit under the roof"""
    cx, cy, cz = P["center"]
    L, W = float(P["length"]), float(P["width"])
    
    # Calculate inset boundaries
    xL, xR = cx - L*0.5, cx + L*0.5
    yB, yT = cy - W*0.5, cy + W*0.5
    xL2, xR2 = xL + P["offset_A"], xR - P["offset_B"]
    yB2, yT2 = yB + P["offset_C"], yT - P["offset_D"]
    
    # Create inset rectangle and extrude
    rectB = add_rect_by_bounds(xL2, xR2, yB2, yT2, cz)
    walls_id = rs.ExtrudeCurveStraight(rectB, (0,0,0), (0,0,float(P["b_extrude"])))
    
    # Split walls with roof and keep bottom piece
    pieces = rs.SplitBrep(walls_id, roof_id, True)
    if not pieces or len(pieces) < 2: 
        print("Warning: Walls did not split with roof")
        return rectB, walls_id
    
    # Find piece with lowest maximum Z (bottom piece)
    walls_bottom = min(pieces, key=zmax_of)
    
    # Delete other pieces
    for pid in pieces:
        if pid != walls_bottom: 
            rs.DeleteObject(pid)
    
    return rectB, walls_bottom

def add_aperture_to_walls(walls_id, P, ap_u, ap_width, ap_height, ap_depth):
    """Cut aperture opening in walls using OFFSET rectangle boundary (where walls actually exist)"""
    cx, cy, cz = P["center"]
    L, W = float(P["length"]), float(P["width"])
    
    # Calculate the SAME offset boundaries used for wall creation
    xL, xR = cx - L*0.5, cx + L*0.5
    yB, yT = cy - W*0.5, cy + W*0.5
    xL2, xR2 = xL + P["offset_A"], xR - P["offset_B"]
    yB2, yT2 = yB + P["offset_C"], yT - P["offset_D"]
    
    # Create OFFSET rectangle boundary for aperture positioning (matches actual wall location)
    offset_boundary_id = add_rect_by_bounds(xL2, xR2, yB2, yT2, cz)
    offset_edges = [e for e in rs.ExplodeCurves(offset_boundary_id, True) or [] if e]
    boundary_ap_list = rs.JoinCurves(offset_edges, True)
    
    if not boundary_ap_list:
        raise RuntimeError("Failed to build offset aperture boundary")
    boundary_ap = boundary_ap_list[0]
    
    # Create aperture cutting box positioned along ACTUAL wall perimeter
    box_id, box_bb = make_aperture_box_on_perimeter(
        boundary_ap, ap_u, ap_width, ap_height, ap_depth
    )
    
    # Calculate bounding box limits
    bb_min_x = min(p.X for p in box_bb)
    bb_max_x = max(p.X for p in box_bb)
    bb_min_y = min(p.Y for p in box_bb)
    bb_max_y = max(p.Y for p in box_bb)
    bb_min_z = min(p.Z for p in box_bb)
    bb_max_z = max(p.Z for p in box_bb)
    
    # Split walls with aperture box
    split_pieces = rs.SplitBrep(walls_id, box_id, True)
    
    # Clean up
    rs.DeleteObject(box_id)
    rs.DeleteObject(boundary_ap)
    rs.DeleteObject(offset_boundary_id)
    for e in offset_edges:
        try:
            if e: rs.DeleteObject(e)
        except: pass
    
    # Remove aperture volume and keep wall pieces
    if split_pieces and len(split_pieces) >= 2:
        survivors = []
        for pid in split_pieces:
            bb = rs.BoundingBox(pid)
            c = bbox_center(bb)
            if (c and (bb_min_x <= c.X <= bb_max_x) and 
                (bb_min_y <= c.Y <= bb_max_y) and 
                (bb_min_z <= c.Z <= bb_max_z)):
                rs.DeleteObject(pid)  # Delete aperture volume
            else:
                survivors.append(pid)  # Keep wall pieces
        
        if survivors:
            final_walls = survivors[0] if len(survivors)==1 else rs.JoinSurfaces(survivors, True)
            return final_walls
    
    return walls_id

def make_skylight_cutting_box(roof_id, P):
    """
    **Create skylight cutting box positioned on roof surface**
    
    Projects XY ground position vertically up to intersect with angled roof,
    then creates CAPPED cutting box at that roof intersection point.
    """
    
    # === STEP 1: CALCULATE GROUND POSITION ===
    cx, cy, cz = P["center"]
    skylight_ground_x = cx + P["skylight_x"]  # X position on ground
    skylight_ground_y = cy + P["skylight_y"]  # Y position on ground
    ground_point = Rhino.Geometry.Point3d(skylight_ground_x, skylight_ground_y, cz)
    
    # === STEP 2: PROJECT UP TO ROOF SURFACE ===
    # Cast ray upward from ground point to find intersection with angled roof
    upward_ray = Rhino.Geometry.Ray3d(ground_point, Rhino.Geometry.Vector3d(0, 0, 1))
    roof_brep = rs.coercebrep(roof_id)
    
    # Find intersection with roof surface
    intersections = Rhino.Geometry.Intersect.Intersection.RayShoot([roof_brep], upward_ray, 1)
    
    if intersections and len(intersections) > 0:
        hit_point = intersections[0]
        if hasattr(hit_point, 'Point'):
            roof_position = hit_point.Point
        else:
            roof_position = hit_point
        print("Skylight positioned on roof at height {:.1f}".format(roof_position.Z))
    else:
        # Fallback: estimate roof height and position skylight there
        estimated_roof_height = cz + P["ridge_height"] + 5.0  # Slightly above estimated roof
        roof_position = Rhino.Geometry.Point3d(skylight_ground_x, skylight_ground_y, estimated_roof_height)
        print("Warning: No roof intersection found - using estimated position at height {:.1f}".format(estimated_roof_height))
    
    # === STEP 3: CREATE CLOSED SKYLIGHT BOX ===
    half_length = P["skylight_length"] * 0.5
    half_width = P["skylight_width"] * 0.5
    
    # Calculate top and bottom Z positions
    top_z = roof_position.Z + 5.0
    bottom_z = roof_position.Z - P["skylight_height"]
    
    # Create skylight corners at TOP of box
    top_corners = [
        (roof_position.X - half_length, roof_position.Y - half_width, top_z),
        (roof_position.X + half_length, roof_position.Y - half_width, top_z),
        (roof_position.X + half_length, roof_position.Y + half_width, top_z),
        (roof_position.X - half_length, roof_position.Y + half_width, top_z),
        (roof_position.X - half_length, roof_position.Y - half_width, top_z)  # Close
    ]
    
    # Create skylight corners at BOTTOM of box
    bottom_corners = [
        (roof_position.X - half_length, roof_position.Y - half_width, bottom_z),
        (roof_position.X + half_length, roof_position.Y - half_width, bottom_z),
        (roof_position.X + half_length, roof_position.Y + half_width, bottom_z),
        (roof_position.X - half_length, roof_position.Y + half_width, bottom_z),
        (roof_position.X - half_length, roof_position.Y - half_width, bottom_z)  # Close
    ]
    
    # Create top and bottom surfaces
    top_curve = rs.AddPolyline(top_corners)
    bottom_curve = rs.AddPolyline(bottom_corners)
    top_surface = rs.AddPlanarSrf(top_curve)[0]
    bottom_surface = rs.AddPlanarSrf(bottom_curve)[0]
    
    # Create side surfaces by lofting between top and bottom
    # Create individual edge curves between corresponding corners
    side_surfaces = []
    for i in range(4):  # 4 sides
        # Get corner indices (wrapping around)
        corner1_top = top_corners[i]
        corner2_top = top_corners[i + 1]
        corner1_bottom = bottom_corners[i]
        corner2_bottom = bottom_corners[i + 1]
        
        # Create rectangular side surface
        side_corners = [corner1_bottom, corner2_bottom, corner2_top, corner1_top, corner1_bottom]
        side_curve = rs.AddPolyline(side_corners)
        side_surface = rs.AddPlanarSrf(side_curve)[0]
        side_surfaces.append(side_surface)
        rs.DeleteObject(side_curve)  # Clean up curve
    
    # Join all surfaces into closed solid
    all_surfaces = [top_surface, bottom_surface] + side_surfaces
    skylight_box = rs.JoinSurfaces(all_surfaces, True)
    
    # Clean up individual surfaces
    rs.DeleteObject(top_curve)
    rs.DeleteObject(bottom_curve)
    for surface in all_surfaces:
        if rs.IsObject(surface):
            rs.DeleteObject(surface)
    
    # Ensure we have a valid solid
    if not skylight_box:
        print("Warning: Failed to create closed skylight box - using fallback method")
        # Fallback to simple extrusion method
        simple_profile = rs.AddPolyline([
            (roof_position.X - half_length, roof_position.Y - half_width, roof_position.Z),
            (roof_position.X + half_length, roof_position.Y - half_width, roof_position.Z),
            (roof_position.X + half_length, roof_position.Y + half_width, roof_position.Z),
            (roof_position.X - half_length, roof_position.Y + half_width, roof_position.Z),
            (roof_position.X - half_length, roof_position.Y - half_width, roof_position.Z)
        ])
        skylight_box = rs.ExtrudeCurveStraight(simple_profile, 
                                               (0, 0, 5), 
                                               (0, 0, -P["skylight_height"]))
        rs.DeleteObject(simple_profile)
    
    print("Created closed skylight cutting box")
    
    # Return box and bounding box for cutting operations
    bb = rs.BoundingBox(skylight_box)
    return skylight_box, bb

def add_skylight_to_building(roof_id, walls_id, P):
    """
    **Cut skylight opening in roof and any overlapping walls**
    
    Uses same proven cutting logic as doors - bounding box intersection method.
    """
    
    # === STEP 1: CREATE SKYLIGHT CUTTING BOX ===
    skylight_box, box_bb = make_skylight_cutting_box(roof_id, P)
    
    print("Cutting skylight at position ({:.1f}, {:.1f})...".format(P["skylight_x"], P["skylight_y"]))
    
    # Calculate bounding box limits (SAME METHOD AS DOORS)
    bb_min_x = min(p.X for p in box_bb)
    bb_max_x = max(p.X for p in box_bb)
    bb_min_y = min(p.Y for p in box_bb)
    bb_max_y = max(p.Y for p in box_bb)
    bb_min_z = min(p.Z for p in box_bb)
    bb_max_z = max(p.Z for p in box_bb)
    
    # === STEP 2: CUT ROOF WITH SKYLIGHT (SAME LOGIC AS DOORS) ===
    roof_pieces = rs.SplitBrep(roof_id, skylight_box, True)
    final_roof = roof_id  # Default to original
    
    if roof_pieces and len(roof_pieces) >= 2:
        print("Roof split into {} pieces".format(len(roof_pieces)))
        roof_survivors = []
        
        # Use SAME piece selection logic as doors
        for pid in roof_pieces:
            bb = rs.BoundingBox(pid)
            c = bbox_center(bb)
            
            # If piece center is within skylight bounds, delete it (it's the opening)
            if (c and (bb_min_x <= c.X <= bb_max_x) and 
                (bb_min_y <= c.Y <= bb_max_y) and 
                (bb_min_z <= c.Z <= bb_max_z)):
                rs.DeleteObject(pid)  # Delete skylight volume from roof
                print("Removed roof piece in skylight area")
            else:
                roof_survivors.append(pid)  # Keep roof pieces
        
        # Join remaining roof pieces (SAME AS DOORS)
        if roof_survivors:
            if len(roof_survivors) == 1:
                final_roof = roof_survivors[0]
            else:
                joined_roof = rs.JoinSurfaces(roof_survivors, True)
                final_roof = joined_roof if joined_roof else roof_survivors[0]
            rs.DeleteObject(roof_id)  # Remove original roof
            print("Skylight cut in roof successful - {} pieces remaining".format(len(roof_survivors)))
        else:
            print("ERROR: No roof pieces survived - keeping original roof")
            final_roof = roof_id  # Keep original if all pieces deleted
    else:
        print("Roof splitting failed - skylight may be outside roof bounds")
    
    # === STEP 3: CUT WALLS WITH SKYLIGHT (MORE CONSERVATIVE FOR LARGE SKYLIGHTS) ===
    wall_pieces = rs.SplitBrep(walls_id, skylight_box, True)
    final_walls = walls_id  # Default to original
    
    if wall_pieces and len(wall_pieces) >= 2:
        print("Walls split into {} pieces".format(len(wall_pieces)))
        wall_survivors = []
        
        # Calculate skylight center for more precise wall cutting
        skylight_center_x = P["center"][0] + P["skylight_x"]
        skylight_center_y = P["center"][1] + P["skylight_y"]
        
        # Use MORE CONSERVATIVE logic for walls (only delete pieces very close to skylight center)
        for pid in wall_pieces:
            bb = rs.BoundingBox(pid)
            c = bbox_center(bb)
            
            # Only delete wall pieces if their center is VERY close to skylight center
            # (much more conservative than roof cutting)
            if c:
                distance_to_skylight_center = ((c.X - skylight_center_x)**2 + (c.Y - skylight_center_y)**2)**0.5
                max_delete_distance = min(P["skylight_length"], P["skylight_width"]) * 0.3  # Only 30% of skylight size
                
                if distance_to_skylight_center < max_delete_distance:
                    rs.DeleteObject(pid)  # Delete only pieces very close to skylight center
                    print("Removed wall piece very close to skylight center")
                else:
                    wall_survivors.append(pid)  # Keep most wall pieces
            else:
                wall_survivors.append(pid)  # Keep if can't determine center
        
        # Join remaining wall pieces (SAME AS DOORS)
        if wall_survivors:
            if len(wall_survivors) == 1:
                final_walls = wall_survivors[0]
            else:
                joined_walls = rs.JoinSurfaces(wall_survivors, True)
                final_walls = joined_walls if joined_walls else wall_survivors[0]
            rs.DeleteObject(walls_id)  # Remove original walls
            print("Skylight cut in walls successful - {} pieces remaining".format(len(wall_survivors)))
        else:
            print("Warning: All wall pieces would be deleted - keeping original walls")
            final_walls = walls_id  # Keep original if all pieces would be deleted
    else:
        print("Walls not affected by skylight")
    
    # === STEP 4: CLEANUP ===
    rs.DeleteObject(skylight_box)
    
    return final_roof, final_walls
    """Cut aperture opening in walls using OFFSET rectangle boundary (where walls actually exist)"""
    cx, cy, cz = P["center"]
    L, W = float(P["length"]), float(P["width"])
    
    # Calculate the SAME offset boundaries used for wall creation
    xL, xR = cx - L*0.5, cx + L*0.5
    yB, yT = cy - W*0.5, cy + W*0.5
    xL2, xR2 = xL + P["offset_A"], xR - P["offset_B"]
    yB2, yT2 = yB + P["offset_C"], yT - P["offset_D"]
    
    # Create OFFSET rectangle boundary for aperture positioning (matches actual wall location)
    offset_boundary_id = add_rect_by_bounds(xL2, xR2, yB2, yT2, cz)
    offset_edges = [e for e in rs.ExplodeCurves(offset_boundary_id, True) or [] if e]
    boundary_ap_list = rs.JoinCurves(offset_edges, True)
    
    if not boundary_ap_list:
        raise RuntimeError("Failed to build offset aperture boundary")
    boundary_ap = boundary_ap_list[0]
    
    # Create aperture cutting box positioned along ACTUAL wall perimeter
    box_id, box_bb = make_aperture_box_on_perimeter(
        boundary_ap, ap_u, ap_width, ap_height, ap_depth
    )
    
    # Calculate bounding box limits
    bb_min_x = min(p.X for p in box_bb)
    bb_max_x = max(p.X for p in box_bb)
    bb_min_y = min(p.Y for p in box_bb)
    bb_max_y = max(p.Y for p in box_bb)
    bb_min_z = min(p.Z for p in box_bb)
    bb_max_z = max(p.Z for p in box_bb)
    
    # Split walls with aperture box
    split_pieces = rs.SplitBrep(walls_id, box_id, True)
    
    # Clean up
    rs.DeleteObject(box_id)
    rs.DeleteObject(boundary_ap)
    rs.DeleteObject(offset_boundary_id)
    for e in offset_edges:
        try:
            if e: rs.DeleteObject(e)
        except: pass
    
    # Remove aperture volume and keep wall pieces
    if split_pieces and len(split_pieces) >= 2:
        survivors = []
        for pid in split_pieces:
            bb = rs.BoundingBox(pid)
            c = bbox_center(bb)
            if (c and (bb_min_x <= c.X <= bb_max_x) and 
                (bb_min_y <= c.Y <= bb_max_y) and 
                (bb_min_z <= c.Z <= bb_max_z)):
                rs.DeleteObject(pid)  # Delete aperture volume
            else:
                survivors.append(pid)  # Keep wall pieces
        
        if survivors:
            final_walls = survivors[0] if len(survivors)==1 else rs.JoinSurfaces(survivors, True)
            return final_walls
    
    return walls_id

# =============================================================================
# MAIN EXECUTION - Complete Workflow Process
# =============================================================================

def build_specs_text(P, run_index):
    # ASCII-only: 'deg', 'x', plain quotes
    return (
        "House #{}\n"
        "Roof: +{:.0f} deg / {:.0f} deg @ z{:.1f}\n"
        "Size: L{:.1f} x W{:.1f}\n"
        "Walls: h{:.1f}\n"
        "Doors: w{:.1f}x{:.1f} @u{:.2f}, w{:.1f}x{:.1f} @u{:.2f}\n"
        "Skylight: L{:.1f}xW{:.1f} @({:.1f},{:.1f},z{:.1f})\n"
        "Additional Surface C: Center({:.1f},{:.1f},{:.1f}) L{:.1f} x W{:.1f}"
    ).format(
        run_index + 1,
        float(P.get("rotate_deg_pos", 0.0)), float(P.get("rotate_deg_neg", 0.0)), float(P.get("ridge_height", 0.0)),
        float(P.get("length", 0.0)), float(P.get("width", 0.0)),
        float(P.get("b_extrude", 0.0)),
        float(P.get("ap_width", 0.0)), float(P.get("ap_height", 0.0)), float(P.get("ap_u", 0.0)),
        float(P.get("ap2_width", 0.0)), float(P.get("ap2_height", 0.0)), float(P.get("ap2_u", 0.0)),
        float(P.get("skylight_length", 0.0)), float(P.get("skylight_width", 0.0)),
        float(P.get("skylight_x", 0.0)), float(P.get("skylight_y", 0.0)), float(P.get("skylight_z", 0.0)),
        float(P.get("C_center", (0.0,0.0,0.0))[0]), float(P.get("C_center", (0.0,0.0,0.0))[1]), float(P.get("C_center", (0.0,0.0,0.0))[2]),
        float(P.get("C_length", 0.0)), float(P.get("C_width", 0.0))
    )

def place_specs_label_dot(obj_ids, text, pad=2.0):
    """Place a TextDot just 'in front' (-Y) of the combined bbox."""
    if not obj_ids: return None
    bbs = [rs.BoundingBox(oid) for oid in obj_ids if rs.IsObject(oid)]
    pts = [p for bb in bbs if bb for p in bb]
    if not pts: return None
    xs = [p.X for p in pts]; ys = [p.Y for p in pts]; zs = [p.Z for p in pts]
    xmin, xmax = min(xs), max(xs)
    ymin = min(ys); zmin = min(zs)
    px = 0.5*(xmin + xmax)
    py = ymin - pad
    pz = zmin + 0.2
    pt = Rhino.Geometry.Point3d(px, py, pz)
    tid = rs.AddTextDot(text, pt)
    if tid: rs.ObjectName(tid, "SpecsLabel")
    return tid

def generate_roof_structure():
    """Main function creating complete building with roof, walls, doors, and skylight"""
    P = PARAMS
    print("Generating building with roof, walls, dual doors, and skylight...")
    
    # === STEP 1: CREATE BASE RECTANGLE ===
    cx, cy, cz = P["center"]
    L, W = float(P["length"]), float(P["width"])
    rectA_id = add_centered_rectangle((cx,cy,cz), L, W)
    
    # === STEP 2: CREATE GABLE ROOF ===
    roof_id = create_gable_roof(rectA_id, P)
    
    # === STEP 3: ELEVATE ROOF TO FINAL HEIGHT ===
    move_xf = Rhino.Geometry.Transform.Translation(0, 0, float(P["ridge_height"]))
    sc.doc.Objects.Transform(roof_id, move_xf, True)
    
    # === STEP 4: CREATE WALLS UNDER ROOF ===
    wall_base, walls_bottom = create_walls_under_roof(P, roof_id)
    
    # === STEP 5: ADD FIRST DOOR ===
    walls_with_door1 = add_aperture_to_walls(
        walls_bottom, P, P["ap_u"], P["ap_width"], P["ap_height"], P["ap_depth"]
    )
    
    # === STEP 6: ADD SECOND DOOR ===
    walls_with_both_doors = add_aperture_to_walls(
        walls_with_door1, P, P["ap2_u"], P["ap2_width"], P["ap2_height"], P["ap2_depth"]
    )
    
    # === STEP 7: ADD SKYLIGHT (CUTS ROOF AND OVERLAPPING WALLS) ===
    final_roof, final_walls = add_skylight_to_building(roof_id, walls_with_both_doors, P)
    
    # === STEP 8: CREATE ADDITIONAL SURFACE ===
    rectC_id = add_centered_rectangle(P["C_center"], P["C_length"], P["C_width"])
    rectC_srf = rs.AddPlanarSrf(rectC_id)[0]
    rs.ObjectName(rectC_srf, "RectC_Surface")
    
    # === STEP 9: CLEANUP ===
    cleanup_objects = [rectA_id, wall_base, rectC_id]
    for obj in cleanup_objects:
        try:
            if obj and rs.IsObject(obj):
                rs.DeleteObject(obj)
        except:
            pass
    
    # === STEP 10: SELECT AND DISPLAY RESULTS ===
    final_objects = [final_roof, final_walls, rectC_srf]
    valid_objects = [obj for obj in final_objects if obj and rs.IsObject(obj)]
    
        # ---- Option B: step each run along +X by (length + GAP) ----
    GAP = 5.0  # spacing between runs (model units)
    run_index = int(sc.sticky.get("gable_run_index", 0))
    dx = run_index * (float(P["length"]) + GAP)

    if valid_objects and dx != 0.0:
        xform = Rhino.Geometry.Transform.Translation(dx, 0.0, 0.0)
        for oid in valid_objects:
            sc.doc.Objects.Transform(oid, xform, True)

    sc.sticky["gable_run_index"] = run_index + 1
    # -------------------------------------------------------------
    

    # Group label with building so they move together if needed
    specs_text = build_specs_text(PARAMS, run_index) # Pass run_index here
    #label_id = place_specs_label_dot(valid_objects, specs_text, pad=2.0)
    # Log the same specs to file
    log_specs_to_file(specs_text, "C:/temp/building_specs_log.txt", auto_open=True)

# Optionally group label with geometry
    if label_id:
        gid = rs.AddGroup("Specs")
        rs.AddObjectsToGroup(valid_objects + [label_id], gid)
        if label_id:
            gname = "BuildingRun_{:03d}".format(run_index)
            gid = rs.AddGroup(gname)
            rs.AddObjectsToGroup(valid_objects + [label_id], gid)
    
    if valid_objects:
        rs.SelectObjects(valid_objects)
    rs.Redraw()
    
    print("Building generation complete!")
    print("- Gable roof: {:.1f} and {:.1f} angles at {:.1f} height".format(
        P["rotate_deg_pos"], P["rotate_deg_neg"], P["ridge_height"]))
    print("- Walls: {:.1f} height with dual doors".format(P["b_extrude"]))
    print("- Door 1: {:.1f}W x {:.1f}H at position {:.2f}".format(
        P["ap_width"], P["ap_height"], P["ap_u"]))
    print("- Door 2: {:.1f}W x {:.1f}H at position {:.2f}".format(
        P["ap2_width"], P["ap2_height"], P["ap2_u"]))
    print("- Skylight: {:.1f}L x {:.1f}W at ({:.1f}, {:.1f})".format(
        P["skylight_length"], P["skylight_width"], P["skylight_x"], P["skylight_y"]))

# =============================================================================
# EXECUTION
# =============================================================================

if __name__ == "__main__":
    generate_roof_structure()
