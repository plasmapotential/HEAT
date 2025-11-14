#meshOpsClass.py
#Description:   HEAT module to create non-VTK and non-STL mesh objects (.gltf, .glb, .usd, etc)
#Engineer:      T Looby
#Date:          20250827

import numpy as np
from pygltflib import (
    GLTF2, Scene, Node, Mesh, Buffer, BufferView, Accessor, Asset, Primitive,
    ELEMENT_ARRAY_BUFFER, ARRAY_BUFFER, Material
)
from struct import pack
import os
from pxr import Usd, UsdGeom, Sdf, Vt, Gf
import plotly.colors as pc

class meshOps:
    def __init__(self):
        return

    def initializeMeshScalar(self, mesh, scalar, label):
        """
        initializes a HEAT mesh object which requires a
        freecad mesh object and a scalar numpy array defined on the mesh

        scalar should be the length of mesh triangles in the mesh

        label is the name of the scalar that will be in the colorbar in paraview
        """
        self.mesh = mesh
        self.scalar = scalar
        self.Npts = len(scalar)
        self.label = label
        return


    def values_to_plotly_rgba(self, values, colorscale="Jet"):
        """
        Map normalized values in [0,1] to RGBA uint8 using a Plotly colorscale.
        Accepts both named Plotly colorscales ("Viridis", "Inferno", etc.)
        and custom scales ([[0,"blue"],[1,"red"]]).
        """

        # get the colorscale
        if isinstance(colorscale, str):
            cs = pc.get_colorscale(colorscale)
        else:
            cs = colorscale

        def parse_plotly_color(c):
            """Convert a Plotly color string ("rgb()", "rgba()", "#hex") to (r,g,b)."""
            if isinstance(c, str):
                if c.startswith("#"):  # hex
                    return pc.hex_to_rgb(c)
                elif c.startswith("rgb"):  # rgb or rgba
                    nums = pc.unlabel_rgb(c)
                    return tuple(int(x) for x in nums[:3])  # drop alpha if present
            elif isinstance(c, (list, tuple)) and len(c) >= 3:
                return tuple(int(x) for x in c[:3])
            raise ValueError(f"Unsupported color format: {c}")

        def interp(val):
            """Interpolate a single value against the colorscale."""
            for i in range(1, len(cs)):
                if val <= cs[i][0]:
                    f0, c0 = cs[i-1]
                    f1, c1 = cs[i]
                    w = (val - f0) / (f1 - f0) if f1 > f0 else 0
                    r0, g0, b0 = parse_plotly_color(c0)
                    r1, g1, b1 = parse_plotly_color(c1)
                    r = int((1-w)*r0 + w*r1)
                    g = int((1-w)*g0 + w*g1)
                    b = int((1-w)*b0 + w*b1)
                    return [r, g, b, 255]
            # fallback: last color
            r, g, b = parse_plotly_color(cs[-1][1])
            return [r, g, b, 255]

        return np.array([interp(v) for v in values], dtype=np.uint8)



    def writeMeshGLB(self, outFile, clim=None, cmap="Jet", embed_scalar_as_custom=True, unlit=False):
        """
        Write a valid .glb with:
          - POSITION (float32)
          - COLOR_0 (RGBA u8, normalized) from heat-flux for out-of-the-box coloring
          - (optional) _HEAT_q (float32 SCALAR) per-vertex to preserve exact MW/m^2

        Requires: pip install pygltflib
        """
        # --------------------------
        # Build per-vertex arrays
        # --------------------------
        verts = []
        indices = []
        facet_vals = np.asarray(self.scalar, dtype=np.float32)
        nF = len(self.mesh.Facets)

        for i, facet in enumerate(self.mesh.Facets):
            base = len(verts)
            for j in range(3):
                x, y, z = facet.Points[j]
                verts.append((float(x), float(y), float(z)))
            # triangle uses the 3 freshly added vertices
            indices.extend([base, base + 1, base + 2])

        V = np.asarray(verts, dtype=np.float32)             # (Nvert, 3)
        I = np.asarray(indices, dtype=np.uint32).reshape(-1) # (Ntri*3,)
        Nvert = V.shape[0]
        assert I.max(initial=0) < Nvert

        # per-vertex heat flux: replicate facet scalar to each of its 3 verts
        HF = np.repeat(facet_vals, 3).astype(np.float32)     # (Nvert,)

        # --------------------------
        # Map HF -> vertex RGBA u8 (Plotly colormap)
        # --------------------------
        if HF.size == 0:
            t = np.zeros(0, dtype=np.float32)
        else:
            if clim is None:
                vmin = float(np.nanmin(HF))
                vmax = float(np.nanmax(HF))
                if not np.isfinite(vmin): vmin = 0.0
                if not np.isfinite(vmax) or vmax <= vmin: vmax = vmin + 1.0
            else:
                vmin, vmax = map(float, clim)
                if vmax <= vmin: vmax = vmin + 1.0
            t = np.clip((HF - vmin) / (vmax - vmin), 0.0, 1.0).astype(np.float32)

        C_u8 = self.values_to_plotly_rgba(t, cmap)

        # (optional) exact floats we may embed as custom attribute
        Q_f32 = HF

        # --------------------------
        # Pack one binary buffer with 4-byte padding between views
        # --------------------------
        def pad4(n): return (4 - (n % 4)) % 4

        idx_blob = I.tobytes()
        idx_pad  = b"\x00" * pad4(len(idx_blob))

        pos_blob = V.astype(np.float32).tobytes()
        pos_pad  = b"\x00" * pad4(len(pos_blob))

        col_blob = C_u8.tobytes()
        col_pad  = b"\x00" * pad4(len(col_blob))

        offset_idx = 0
        offset_pos = offset_idx + len(idx_blob) + len(idx_pad)
        offset_col = offset_pos + len(pos_blob) + len(pos_pad)

        full_blob = idx_blob + idx_pad + pos_blob + pos_pad + col_blob + col_pad

        # optional custom attribute blob
        if embed_scalar_as_custom and Q_f32.size:
            q_blob = Q_f32.tobytes()
            q_pad  = b"\x00" * pad4(len(q_blob))
            offset_q = len(full_blob)
            full_blob += q_blob + q_pad
        else:
            offset_q = None

        # --------------------------
        # Build glTF structure
        # --------------------------

        # BufferViews
        bv_idx = BufferView(buffer=0, byteOffset=offset_idx, byteLength=len(idx_blob), target=ELEMENT_ARRAY_BUFFER)
        bv_pos = BufferView(buffer=0, byteOffset=offset_pos, byteLength=len(pos_blob), target=ARRAY_BUFFER)
        bv_col = BufferView(buffer=0, byteOffset=offset_col, byteLength=len(col_blob), target=ARRAY_BUFFER)
        bufferViews = [bv_idx, bv_pos, bv_col]
        if offset_q is not None:
            bv_q = BufferView(buffer=0, byteOffset=offset_q, byteLength=len(q_blob), target=ARRAY_BUFFER)
            bufferViews.append(bv_q)

        # Accessors
        acc_idx = Accessor(bufferView=0, byteOffset=0, componentType=5125,  # UNSIGNED_INT
                           count=I.size, type="SCALAR",
                           max=[int(I.max())] if I.size else [0],
                           min=[int(I.min())] if I.size else [0])
        acc_pos = Accessor(bufferView=1, byteOffset=0, componentType=5126,  # FLOAT
                           count=Nvert, type="VEC3",
                           max=list(map(float, V.max(axis=0))) if Nvert else [0,0,0],
                           min=list(map(float, V.min(axis=0))) if Nvert else [0,0,0])
        acc_col = Accessor(bufferView=2, byteOffset=0, componentType=5121,  # UNSIGNED_BYTE
                           count=Nvert, type="VEC4", normalized=True)

        accessors = [acc_idx, acc_pos, acc_col]
        if offset_q is not None:
            acc_q = Accessor(bufferView=3, byteOffset=0, componentType=5126,  # FLOAT
                             count=Nvert, type="SCALAR",
                             max=[float(np.nanmax(Q_f32))], min=[float(np.nanmin(Q_f32))])
            accessors.append(acc_q)

        # Primitive attributes
        attrs = {"POSITION": 1, "COLOR_0": 2}
        if offset_q is not None:
            attrs[self.label] = 3

        # --- Material handling ---
        if unlit:
            mat = Material(
                pbrMetallicRoughness={}, 
                extensions={"KHR_materials_unlit": {}}
            )
            extensionsUsed = ["KHR_materials_unlit"]
        else:
            mat = Material(pbrMetallicRoughness={})
            extensionsUsed = []

        prim = Primitive(
            attributes=attrs,
            indices=0,
            mode=4,  # TRIANGLES
            material=0
        )

        mesh = Mesh(primitives=[prim])
        node = Node(mesh=0, name="HEATMesh")
        scene = Scene(nodes=[0])

        gltf = GLTF2(
            asset=Asset(version="2.0"),
            buffers=[Buffer(byteLength=len(full_blob))],
            bufferViews=bufferViews,
            accessors=accessors,
            meshes=[mesh],
            nodes=[node],
            scenes=[scene],
            scene=0,
            materials=[mat],
            extensionsUsed=extensionsUsed
        )

        # Embed the binary and save
        gltf.set_binary_blob(full_blob)

        # ensure .glb
        base, ext = os.path.splitext(outFile)
        if ext.lower() != ".glb":
            outFile = base + ".glb"

        gltf.save_binary(outFile)
        return outFile



    def writeMeshUSD(self,
                     outFile,
                     clim=None,
                     colormap='grayscale',
                     up_axis='Z',            # 'Z' for tokamak coords, use 'Y' if you prefer Omniverse default
                     meters_per_unit=1.0,    # set stage scale (1.0 = meters)
                     double_sided=True,
                     write_normals=True):
        """
        Write a USD file with:
          - topology (triangles),
          - primvars:displayColor (vertex colors from HF for instant visualization),
          - primvars:heatFlux (vertex float array with exact MW/m^2 values).

        Requires: pip install usd-core  (gives pxr.*) or a Python that ships with USD.
        """
        # ----------------------------
        # Build flattened vertices/indices and per-vertex scalars
        # ----------------------------
        pts = []
        fvi = []
        fvc = []
        scalars = []

        for i, facet in enumerate(self.mesh.Facets):
            base = len(pts)
            for j in range(3):
                x, y, z = facet.Points[j]
                pts.append((float(x), float(y), float(z)))
                scalars.append(float(self.scalar[i]))  # replicate facet scalar to the 3 vertices
            fvi.extend([base, base + 1, base + 2])
            fvc.append(3)

        if not pts:
            raise ValueError("No facets/points found.")

        # numpy views for convenience
        P = np.asarray(pts, dtype=np.float32)
        HF = np.asarray(scalars, dtype=np.float32)

        # ----------------------------
        # Colormap for displayColor (optional but recommended)
        # ----------------------------
        if clim is None:
            vmin = float(np.nanmin(HF)) if HF.size else 0.0
            vmax = float(np.nanmax(HF)) if HF.size else 1.0
            if vmax <= vmin:
                vmax = vmin + 1e-6
        else:
            vmin, vmax = map(float, clim)
            if vmax <= vmin:
                vmax = vmin + 1e-6

        t = np.clip((HF - vmin) / (vmax - vmin), 0.0, 1.0)
        if colormap == 'grayscale':
            C = np.stack([t, t, t], axis=1).astype(np.float32)   # Color3f in 0..1
        else:
            # simple fallback (plug in your LUT here if needed)
            C = np.stack([t, t, t], axis=1).astype(np.float32)

        # ----------------------------
        # Optional normals (vertex; here we do per-face normals since verts are duplicated)
        # ----------------------------
        N = None
        if write_normals:
            N = np.zeros_like(P)
            tris = np.asarray(fvi, dtype=np.int32).reshape(-1, 3)
            a = P[tris[:, 0]]
            b = P[tris[:, 1]]
            c = P[tris[:, 2]]
            n = np.cross(b - a, c - a)
            lens = np.linalg.norm(n, axis=1)
            mask = lens > 0.0
            n[mask] /= lens[mask, None]
            # assign same face normal to each of its 3 vertices
            N[tris[:, 0]] = n
            N[tris[:, 1]] = n
            N[tris[:, 2]] = n

        # ----------------------------
        # Create USD stage & prims
        # ----------------------------
        if not outFile.lower().endswith(('.usd', '.usdc', '.usda')):
            outFile = outFile + '.usdc'

        stage = Usd.Stage.CreateNew(outFile)
        # Stage metadata
        UsdGeom.SetStageMetersPerUnit(stage, float(meters_per_unit))
        UsdGeom.SetStageUpAxis(stage, UsdGeom.Tokens.z if up_axis.upper().startswith('Z') else UsdGeom.Tokens.y)

        # Optional root xform
        UsdGeom.Xform.Define(stage, Sdf.Path('/World'))

        mesh = UsdGeom.Mesh.Define(stage, Sdf.Path('/World/HEATMesh'))
        mesh.CreateSubdivisionSchemeAttr().Set(UsdGeom.Tokens.none)
        mesh.CreateOrientationAttr().Set(UsdGeom.Tokens.rightHanded)
        mesh.CreateDoubleSidedAttr().Set(bool(double_sided))

        # Topology & points
        mesh.CreatePointsAttr(Vt.Vec3fArray([Gf.Vec3f(*p) for p in P.tolist()]))
        mesh.CreateFaceVertexCountsAttr(Vt.IntArray(fvc))
        mesh.CreateFaceVertexIndicesAttr(Vt.IntArray(fvi))

        # Extent (AABB)
        mn = P.min(axis=0)
        mx = P.max(axis=0)
        mesh.CreateExtentAttr(Vt.Vec3fArray([Gf.Vec3f(*mn), Gf.Vec3f(*mx)]))

        # displayColor (what most DCCs/viewers show by default)
        pv_api = UsdGeom.PrimvarsAPI(mesh)
        color_pv = pv_api.CreatePrimvar('displayColor',
                                        Sdf.ValueTypeNames.Color3fArray,
                                        UsdGeom.Tokens.vertex)
        color_pv.Set(Vt.Vec3fArray([Gf.Vec3f(*rgb) for rgb in C.tolist()]))

        # exact heat-flux values as vertex primvar (MW/m^2)
        hf_pv = pv_api.CreatePrimvar('heatFlux',
                                     Sdf.ValueTypeNames.FloatArray,
                                     UsdGeom.Tokens.vertex)
        hf_pv.Set(Vt.FloatArray([float(x) for x in HF.tolist()]))
        # advertise units via custom metadata
        hf_pv.GetAttr().SetMetadata('documentation', self.label)

        # Normals (optional)
        if N is not None:
            n_pv = pv_api.CreatePrimvar('normals',
                                        Sdf.ValueTypeNames.Normal3fArray,
                                        UsdGeom.Tokens.vertex)
            n_pv.Set(Vt.Vec3fArray([Gf.Vec3f(*vv) for vv in N.tolist()]))

        # Save
        stage.GetRootLayer().Save()
        return outFile
