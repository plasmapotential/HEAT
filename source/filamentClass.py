#filamentClass.py
#Description:   HEAT filament module
#Engineer:      T Looby
#Date:          20230228
"""
Class for filament tracing.  Used for ELMs and other filamentary structures
"""
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import scipy.interpolate as interp
import pandas as pd

#HEAT classes
import toolsClass
tools = toolsClass.tools()
import ioClass
IO = ioClass.IO_HEAT()


class filament:

    def __init__(self, rootDir, dataPath, chmod=0o774, UID=-1, GID=-1):
        """
        rootDir is root location of python modules (where dashGUI.py lives)
        dataPath is the location where we write all output to
        """
        self.rootDir = rootDir
        tools.rootDir = self.rootDir
        self.dataPath = dataPath
        tools.dataPath = self.dataPath
        self.chmod = chmod
        self.GID = GID
        self.UID = UID
        return

    def setupNumberFormats(self, tsSigFigs=6, shotSigFigs=6):
        """
        sets up pythonic string number formats for shot and timesteps
        """
        self.tsFmt = "{:."+"{:d}".format(tsSigFigs)+"f}"
        self.shotFmt = "{:0"+"{:d}".format(shotSigFigs)+"d}"
        return

    def readFilamentFile(self, path):
        """
        reads a filament csv input file

        generates a data pd array
        """

        filFile = path + 'filaments.csv'
        self.filData = pd.read_csv(filFile, sep=',', comment='#', skipinitialspace=True)
        return

    def setupFilamentTime(self):
        """
        sets up filament timesteps using data from filament file
        """
        ts = np.array([])
        for row in self.filData:
            tMin = row['tMin[s]']
            tMax = row['tMax[s]']
            dt = row['dt']
            N_dt = int( (tMax-tMin)/dt )
            tNew = np.linspace(tMin,tMax,N_dt+1)
            ts = np.append(ts, tNew)
        
        self.tsFil = ts
        return

    def allowed_class_vars(self):
        """
        Writes a list of recognized class variables to HEAT object
        Used for error checking input files and for initialization

        Here is a list of variables with description:
        testvar         dummy for testing

        """


        self.allowed_vars = [
                            'testVar',
                            ]
        return

    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file
        """

        integers = [                  
                    ]

        floats = [
                  'testVar',
                 ]


        for var in integers:
            if (getattr(self, var) is not None) and (~np.isnan(float(getattr(self, var)))):
                try:
                    setattr(self, var, tools.makeInt(getattr(self, var)))
                except:
                    print("Error with input file var "+var+".  Perhaps you have invalid input values?")
                    log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")
        for var in floats:
            if var is not None:
                if (getattr(self, var) is not None) and (~np.isnan(float(getattr(self, var)))):
                    try:
                        setattr(self, var, tools.makeFloat(getattr(self, var)))
                    except:
                        print("Error with input file var "+var+".  Perhaps you have invalid input values?")
                        log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")

        return
    

    def initializeFilament(self, 
                           rCtr: float, 
                           zCtr: float, 
                           phi: float,
                           sigma_b: float, 
                           sigma_r: float, 
                           sigma_p: float, 
                           E0: float, 
                           ep: object):
        """
        initializes a filament object that can be traced in HEAT

        rCtr: r coordinate center of filament [m]
        zCtr: z coordinate center of filament [m]
        phi: toroidal center of filament [radians]
        sigma_b: characteristic width in parallel direction [m]
        sigma_r: characteristic width in radial (psi) direction [m]
        sigma_p: characteristic width in poloidal direction [m]
                   (note that this is NOT the diamagnetic direction exactly)
        E0: total filament energy at t=0 [MJ]
        ep: an equilParams object

        """

        self.rCtr = rCtr
        self.zCtr = zCtr
        self.phi = phi
        self.sig_b = sigma_b
        self.sig_r = sigma_r
        self.sig_p = sigma_p
        self.E0 = E0
        self.ep = ep
        self.xCtr = rCtr * np.cos(phi)
        self.yCtr = rCtr * np.sin(phi)

        return


    def fluxSurfNorms(self, ep: object, R:np.ndarray, Z:np.ndarray):
        """
        Calculates vectors normal to poloidal flux surfaces at R,Z coordinates

        returns normal vector coordinates in (R, Phi, Z)

        phiNorm is returned as 0
        """       
        if np.isscalar(R):                   
            R = np.array([R])
        if np.isscalar(Z):                   
            Z = np.array([Z])

        Br = ep.BRFunc.ev(R,Z)
        Bz = ep.BZFunc.ev(R,Z)
        Bt = ep.BtFunc.ev(R,Z)
        Bp = np.sqrt(Br**2+Bz**2)

        Bmag = np.sqrt(Bt**2 + Bp**2)


        bt = np.zeros((len(Bt), 3)) #R,t,Z
        bp = np.zeros((len(Bp), 3)) #R,t,Z
    
        #bt[:,1] = Bt / Bmag
        #bp[:,0] = Br / Bmag
        #bp[:,2] = Bz / Bmag

        bt[:,1] = -1.0
        bp[:,0] = Br / Bp
        bp[:,2] = Bz / Bp


        norms = np.cross(bt,bp)


        mag = np.linalg.norm(norms, axis=1)
        r_hat = np.zeros((len(norms), 3))
        r_hat = norms / mag[:, None]

        return r_hat

    def poloidalVectors(self, ep: object, R: np.ndarray ,Z: np.ndarray):
        """
        Calculates vectors tangent to poloidal flux surfaces in RZ plane at R,Z coordinates
        
        returns normal vector coordinates in (R, Phi, Z)

        phiNorm is returned as 0
        """

        if np.isscalar(R):                   
            R = np.array([R])
        if np.isscalar(Z):                   
            Z = np.array([Z])


        Br = ep.BRFunc.ev(R,Z)
        Bz = ep.BZFunc.ev(R,Z)
        Bt = ep.BtFunc.ev(R,Z)
        Bp = np.sqrt(Br**2+Bz**2)

        Bmag = np.sqrt(Bt**2 + Bp**2)

        bp = np.zeros((len(Bp), 3)) #R,phi,Z

        bp[:,0] = Br / Bp
        bp[:,2] = Bz / Bp

        return bp



    def getTraceSection(self, low:float, high:float, trace:np.ndarray):
        """
        returns a section of a trace between low and high
        """
        dist = self.distance(trace)
        idxDist = np.where(np.logical_and(dist >= low, dist <= high))[0]
        return trace[idxDist]
    

    def interpolateTrace(self, traceData:np.ndarray, N:int, addRawData=False):
        """
        given a trace of coordinates, interpolates to N discrete points
        """
        #calculate distance along wall
        dist = self.distance(traceData)
        
        # Resolution we need given the inputs
        resolution = (dist[-1] - dist[0])/float(N)
        print("Resolution: {:f} m".format(resolution))
        # Calculate how many total grid points we need based upon resolution and
        # total distance around curve/wall
        numpoints = int((dist[-1] - dist[0])/resolution)
        # Spline Interpolation (linear) - Make higher resolution wall.
        interpolator = interp.interp1d(dist, traceData, kind='slinear', axis=0)
        alpha = np.linspace(dist[0], dist[-1], numpoints)
        #add in raw data points so that we preserve the corners
        if addRawData == True:
            alpha = np.sort(np.hstack([alpha, dist]))
        
        interpolated_points = interpolator(alpha)
        arr, idx = np.unique(interpolated_points, axis=0, return_index=True)
        interpolated_points = arr[np.argsort(idx)]
     
        return interpolated_points, alpha
    

    def distance(self, traceData:np.ndarray):
        """
        Calculate distance along curve/wall (also called S) of ND points:
        """
        distance = np.cumsum(np.sqrt(np.sum(np.diff(traceData,axis=0)**2,axis=1)))
        distance = np.insert(distance, 0, 0)
        return distance




    def gridPsiThetaDistAtCtr(self, rCtr:float, zCtr:float, multR=10.0, multZ = 10.0):
        """
        generates an RZ grid, and then calculates distances [m] on that grid
        in the psi and theta directions from (self.rCtr, self.zCtr)

        these distances can then be used by gaussian

        multR is multiplier to create R width of grid: R width = multR*sigma_r
        multZ is multiplier to create Z height of grid: Z height = multZ*sigma_z
        """
        ep = self.ep
        sigma_r = self.sig_r
        sigma_p = self.sig_p

        #build rectilinear coordinates around the ctr
        r = np.linspace(rCtr - multR*sigma_r, rCtr + multR*sigma_r, 100)
        z = np.linspace(zCtr - multZ*sigma_p, zCtr + multZ*sigma_p, 100)
        R,Z = np.meshgrid(r,z)

        psiCtr, distPsi, thetaCtr, distTheta = self.fluxCoordDistance(rCtr,zCtr,R,Z)

        #set class variables
        self.distPsi = distPsi
        self.distTheta = distTheta
        self.rGrid = r
        self.zGrid = z
        self.psiCtr = psiCtr
        self.thetaCtr = thetaCtr
        self.gridShape = R.shape
        
        return

    def gaussianAtPts(self, ctrPts: np.ndarray, xyzPts: np.ndarray, t: float, v_r:float):
        """
        calculates a gaussian, centered at ctrPts, evaluated at xyzPts, at time t
        v_r is radial velocity

        saves gaussian values on 3D grid, as well as distances (psi, theta) used to
        evaluate the gaussian        
        """
        ep = self.ep
        sigma_r = self.sig_r
        sigma_p = self.sig_p


        gaussian = np.zeros((xyzPts.shape[:-1]))
        dPsi = np.zeros((xyzPts.shape[:-1]))
        dTheta = np.zeros((xyzPts.shape[:-1]))

        for i,ctr in enumerate(ctrPts):
            rCtr, zCtr, phiCtr = tools.xyz2cyl(ctr[0],ctr[1],ctr[2])
            xyz = xyzPts[i]
            R,Z,phi = tools.xyz2cyl(xyz[:,:,0],xyz[:,:,1],xyz[:,:,2])
            psiCtr, distPsi, thetaCtr, distTheta = self.fluxCoordDistance(rCtr,zCtr,R,Z)
            
            #2D gaussian
            #g = self.gaussian2D(distPsi,distTheta,self.sig_r,self.sig_p,t,v_r,self.E0)
            #3D gaussian
            g = self.gaussian3D(distPsi,distTheta,self.distB[i],self.sig_r,self.sig_p,self.sig_b,t,v_r,self.E0)
            #reshape to meshgrid shape
            gaussian[i,:,:] = g.reshape((len(distPsi), len(distTheta)))
            dPsi[i,:,:] = distPsi.reshape((len(distPsi), len(distTheta)))
            dTheta[i,:,:] = distTheta.reshape((len(distPsi), len(distTheta)))

        self.g_pts = gaussian
        self.distPsi = dPsi
        self.distTheta = dTheta
        return


    def fluxCoordDistance(self, r0:float, z0:float, R:np.ndarray, Z:np.ndarray):
        """
        calculates euclidean distance along flux coordinate surfaces (psi and poloidal)
        from rCtr,zCtr on an R,Z grid

        rCtr and zCtr:  coordinate to calculate flux coordinates at / distance from
        R,Z: meshgrid around rCtr,zCtr

        returns psiCtr, distPsi, thetaCtr, distTheta
        """
        ep = self.ep

        #flux coordinates of ctr
        psiCtr = ep.psiFunc(r0, z0)
        thetaCtr = self.thetaFromRZ(ep, r0, z0)

        #==Find radial coordinates
        #Find coordinate transformation at ctr
        psiaxis = ep.g['psiAxis']
        psiedge = ep.g['psiSep']
        deltaPsi = np.abs(psiedge - psiaxis)
        Bp = ep.BpFunc(r0,z0)
        gradPsi = Bp*r0
        xfm = gradPsi / deltaPsi
        #convert distance in psi to euclidean coordinates
        psiRZ = ep.psiFunc.ev(R,Z)
        distPsi = (psiRZ - psiCtr) / xfm #on grid


        #===Find poloidal coordinates
        thetaRZ = self.thetaFromRZ(ep, R.flatten(), Z.flatten())
        #surfR, surfZ = ep.flux_surface(psiCtr, Npts, thetaRZ[sortIdx])
        surfR = ep.g['lcfs'][:,0]
        surfZ = ep.g['lcfs'][:,1]
        thetaSurf = self.thetaFromRZ(ep, surfR, surfZ)
        #calculate distance along flux surface
        surface = np.vstack((surfR, surfZ)).T
        dSurf = self.distance(surface)
        #interpolator that maps theta back to distance along the flux surface
        interpolator = interp.interp1d(thetaSurf, dSurf, kind='slinear', axis=0)
        dCtr = interpolator(thetaCtr)
        distTheta = interpolator(thetaRZ) - dCtr
        distTheta = distTheta.reshape(R.shape) #on grid

        return psiCtr, distPsi, thetaCtr, distTheta


    def filamentGaussian2D(self, t:float, v_r:float, distPsi:np.ndarray, distTheta:np.ndarray):
        """
        calculates gaussian at time t with radial (psi) velocity v_r       
        """
        #calculate gaussian
        g = self.gaussian2D(distPsi,distTheta,self.sig_r,self.sig_p,t,v_r,self.E0)
        #reshape to meshgrid shape
        g.reshape(distPsi.shape) #on grid
        return g
    

    def gaussian2D(self, dx:np.ndarray, dy:np.ndarray , sigX:float ,sigY:float ,t: float, v_r:float, A=1.0):
        """
        calculates gaussian function with two spatial dimensions and 1 time dimension
        """
        #placeholder for now.  decay time of 500us
        if t==0:
            A=1
        else:
            A = 1 - t/500e-6

        exp1 = np.exp( -1.0*(dx - t*v_r)**2 / (2*sigX**2) )
        exp2 = np.exp( -1.0*(dy**2) / (2*sigY**2) )
        g = A  * exp1 * exp2
        return g

    def filamentGaussian3D(self, t:float, v_r:float, distPsi:np.ndarray, distTheta:np.ndarray, distB:np.ndarray):
        """
        calculates gaussian at time t with radial (psi) velocity v_r in three dimensions       
        """
        #calculate gaussian
        g = self.gaussian3D(distPsi,distTheta,distB,self.sig_r,self.sig_p,self.sig_b,t,v_r,self.E0)
        #reshape to meshgrid shape
        g.reshape(distPsi.shape) #on grid
        return g

    def gaussian3D(self, dx:np.ndarray, dy:np.ndarray, dz:np.ndarray, 
                   sigX:float, sigY:float, sigZ:float, 
                   t:float , v_r:float, A=1.0):
        """
        calculates gaussian function with three spatial dimensions and 1 time dimension
        """
        #placeholder for now.  decay time of 500us
        if t==0:
            A=1
        else:
            A = 1 - t/500e-6

        exp1 = np.exp( -1.0*(dx - t*v_r)**2 / (2*sigX**2) )
        exp2 = np.exp( -1.0*(dy**2) / (2*sigY**2) )
        exp3 = np.exp( -1.0*(dz**2) / (2*sigZ**2) )
        g = A  * exp1 * exp2 * exp3
        return g


    def thetaFromRZ(self, ep:object, R:np.ndarray, Z:np.ndarray):
        """
        calculates poloidal coordinate, theta, at R,Z given an equilibrium (ep) object

        returns theta
        """
        r_hat = self.fluxSurfNorms(ep, R, Z)
        r_hat = np.delete(r_hat, 1, 1)  # delete second column (phi)
        R_hat = np.array([1.0, 0.0])
        theta = np.arccos(np.dot(r_hat, R_hat))
        zMid = ep.g['ZmAxis']
        idx = np.where(Z < zMid)[0]
        theta[idx] *= -1.0
        return theta


    def discretizeFilament(self, N_r: int, N_p: int, N_b: int, Btrace: np.ndarray, 
                           N_sig_r: int, N_sig_p: int, N_sig_b: int):
        """
        takes user input parameters and builds a filament 
        
        N_r: number of discrete points in radial (psi) direction
        N_p: number of discrete points in poloidal direction
        N_b: number of discrete points along B field line
        N_sig_r : # of self.sig_r widths for points to cover in +/- radial direction
        N_sig_p : # of self.sig_p widths for points to cover in +/- poloidal direction
        N_sig_0 : # of self.sig_0 widths for points to cover in +/- parallel direction
        """
        d_b = self.distance(Btrace)
        #find index of filament center
        ctr = np.array([self.xCtr, self.yCtr, self.zCtr])
        idx = np.argmin(np.sum(np.abs(Btrace-ctr), axis=1))      
        d_b = d_b - d_b[idx]
        use = np.where(np.abs(d_b) < N_sig_b*self.sig_b)[0]

        #calculate center coordinates along field line
        self.ctrPts, alpha = self.interpolateTrace(Btrace[use], N_b, addRawData=False)
        self.distB = alpha - np.max(alpha) / 2.0

        #calculate points in parallel (B) direction
        bMax = N_sig_b * self.sig_b
        b_pts = np.linspace(-bMax, bMax, N_b)

        #calculate points in radial (psi) direction
        rMax = N_sig_r * self.sig_r
        r_pts = np.linspace(-rMax, rMax, N_r)

        #calculate points in poloidal direction
        pMax = N_sig_p * self.sig_p
        p_pts = np.linspace(-pMax, pMax, N_p)

        #get radial (psi) vectors for each ctrPt
        R = np.sqrt(self.ctrPts[:,0]**2+self.ctrPts[:,1]**2)
        Z = self.ctrPts[:,2]
        phi = np.arctan2(self.ctrPts[:,1], self.ctrPts[:,0])
        rVec = self.fluxSurfNorms(self.ep, R, Z)

        #convert from R,phi,Z to xyz
        rX, rY, rZ = tools.cyl2xyz(rVec[:,0], rVec[:,2], phi, degrees=False)
        rVecXYZ = np.vstack((rX,rY,rZ)).T
        rMag = np.linalg.norm(rVecXYZ, axis=1)
        rN = rVecXYZ / rMag[:,None]

        #get poloidal vectors for each ctrPt
        pVec = self.poloidalVectors(self.ep, R, Z)
        #convert from R,phi,Z to xyz
        pX, pY, pZ = tools.cyl2xyz(pVec[:,0], pVec[:,2], phi, degrees=False)
        pVecXYZ = np.vstack((pX,pY,pZ)).T
        pMag = np.linalg.norm(pVecXYZ, axis=1)
        pN = pVecXYZ / pMag[:,None]


        #create xyz coordinates for each filament source pt
        xyzPts = np.ones((N_b, N_r, N_p,3))*np.nan
        for i in range(N_b):
            for j in range(N_r):
                for k in range(N_p):
                    xyzPts[i,j,k,:] = self.ctrPts[i] + r_pts[j]*rN[i] - p_pts[k]*pN[i]
        
        self.xyzPts = xyzPts
        self.rN = rN
        self.pN = pN

        return


    def vectorGlyphs(self, ctrs:np.ndarray, vecs:np.ndarray, label:str, path:str, tag:str = None):
        """
        creates a VTP glyph vector for visualization in paraview
        """
        prefix = label
        IO.writeGlyphVTP(ctrs,vecs,label,prefix,path,tag)
        return



    def plotly2DContour(self, x:np.ndarray, y:np.ndarray, z:np.ndarray ,
                        fig:object = None, cs='plasma',zmin=0.0, zmax=1.0, mode='gauss'):
        """
        plots a plotly 2D contour plot 

        z should be in the shape of a np meshgrid whose dimensions match x and y
        """
        if fig == None:
            fig = go.Figure()

        if mode=='psi':
            #plot psi
            fig = go.Figure(data =
                go.Contour(
                    z=z,
                    x=x, # horizontal axis
                    y=y, # vertical axis
                    colorscale=cs,
                    contours_coloring='heatmap',
                    name='psi',
                    showscale=False,
                    ncontours=20,
                    )
             )  
        elif mode=='theta':
            #plot theta
            fig = go.Figure(data =
                go.Contour(
                    z=z,
                    x=x, # horizontal axis
                    y=y, # vertical axis
                    colorscale=cs,
                    contours_coloring='heatmap',
                    name='theta',
                    showscale=False,
                    ncontours=20,
                    )
             )    
        elif mode=='psiDist':
        #plot psiDist
            fig = go.Figure(data =
                go.Contour(
                    z=z,
                    x=y, # horizontal axis
                    y=z, # vertical axis
                    colorscale=cs,
                    contours_coloring='heatmap',
                    name='psiDistance',
                    showscale=False,
                    ncontours=20,
                    )
             )  

        elif mode=='thetaDist':
        #plot thetaDist
            fig = go.Figure(data =
                go.Contour(
                    z=z,
                    x=y, # horizontal axis
                    y=x, # vertical axis
                    colorscale=cs,
                    contours_coloring='heatmap',
                    name='thetaDistance',
                    showscale=False,
                    ncontours=20,
                    )
             )  
        else:
        #plot gaussian
            fig.add_trace(
                go.Contour(
                    z=z,
                    x=x, # horizontal axis
                    y=y, # vertical axis
                    colorscale=cs,
                    contours_coloring='heatmap',
                    name='guassian',
                    showscale=False,
                    ncontours=20,
                    line = dict(width = 0),
                    zmin=zmin,
                    zmax=zmax,
                    )
             )       


        #update the axes
        fig.update_xaxes(range=[min(x), max(x)], title="R[m]")
        fig.update_yaxes(range=[min(y), max(y)], title="Z[m]")
        fig.update_yaxes(scaleanchor = "x",scaleratio = 1,)

        return fig



    def plotlyAddTrace(self, fig:object, trace:np.ndarray, 
                       name:str=None, color:str=None, mode:str=None, markDict:str=None):
        """
        adds a trace to an existing plotly figure

        fig is an existing figure


        trace is rz coordinates of new trace
        """

        if name is None:
            name = 'trace'
        
        if color is None:
            N = len(px.colors.qualitative.Plotly)
            idx = np.random.randint(0,N)
            color = px.colors.qualitative.Plotly[idx]
        
        if mode is None:
            mode="lines"
        
        if mode == "markers" and markDict == None:
            markDict={'symbol':"circle", 'size':16, 'color':'white'}

        fig.add_trace(
            go.Scatter(
                x=trace[:,0],
                y=trace[:,1],
                mode=mode,
                marker=markDict,
                name=name,
                line=dict(
                    color=color
                        ),
                )
                )


        return fig

    def plotlyEQ(self, ep):
        """
        returns a DASH object for use directly in dash app
        """

        psi = ep.g['psiRZ']

        #plot
        levels = sorted(np.append([0.0,0.05,0.1,0.25,0.5,0.75,0.95, 1.0], np.linspace(0.99,psi.max(),15)))
        #psi data
        fig = go.Figure(data =
            go.Contour(
                z=psi,
                x=ep.g['R'], # horizontal axis
                y=ep.g['Z'], # vertical axis
                colorscale='cividis',
                contours_coloring='heatmap',
                name='psi',
                showscale=False,
                ncontours=20,
                )
         )
        #wall in green
        fig.add_trace(
            go.Scatter(
                x=ep.g['wall'][:,0],
                y=ep.g['wall'][:,1],
                mode="markers+lines",
                name="Wall",
                line=dict(
                    color="#19fa1d"
                        ),
                )
                )

        #lcfs from geqdsk
        fig.add_trace(
            go.Scatter(
                x=ep.g['lcfs'][:,0],
                y=ep.g['lcfs'][:,1],
                mode="markers+lines",
                name="Wall",
                line=dict(
                    color="red"
                        ),
                )
                )


        fig.update_layout(
        #    title="{:f} [s]".format(t),
            xaxis_title="R [m]",
            yaxis_title="Z [m]",
            autosize=True,
            #for aspect ratio
            #autosize=False,
            #width=width*1.1,
            #height=height,
            showlegend=False,
            font=dict(
                size=18,
                ),
            margin=dict(
                l=10,
                r=10,
                b=10,
                t=100,
                pad=4
                )
            )

        fig.update_yaxes(scaleanchor = "x",scaleratio = 1,)
        return fig
    
