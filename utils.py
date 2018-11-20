from mayavi import mlab
import numpy as np 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def autoscale_equal_axis(ax):
    [xmin,xmax]=ax.get_xlim()
    [ymin,ymax]=ax.get_ylim()
    [zmin,zmax]=ax.get_zlim()
    axmin = np.min([xmin,ymin,zmin])
    axmax = np.max([xmax,ymax,zmax])
    ax.set_xlim([axmin,axmax])
    ax.set_ylim([axmin,axmax])
    ax.set_zlim([axmin,axmax])
    return 
def plot_3d_positions(positions, radii,ax=None,show_plot=False,with_mayavi=False):

        if ax is None and with_mayavi == False:
            fig = plt.figure(figsize=(16,16))
            ax = fig.add_subplot(111, projection='3d')
           
        #polar coords
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        
        axmin,axmax = +np.inf,-np.inf

        #draw every particle as a surface
        for i in range(len(radii)):
            r = radii[i]
            cx,cy,cz = positions[i,0],positions[i,1],positions[i,2]
            x = r * np.outer(np.cos(u), np.sin(v)) + cx
            y = r * np.outer(np.sin(u), np.sin(v)) + cy
            z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + cz
            for p in [x,y,z]:
                axmin = min(axmin,np.min(p))
                axmax = max(axmax,np.max(p))
            if with_mayavi == False:
                ax.plot_surface(x, y, z,color="red")
                ax.set_xlabel("X [nm]")
                ax.set_ylabel("Y [nm]")
                ax.set_zlabel("Z [nm]")
            elif with_mayavi == True:
                mlab.mesh(x, y, z, representation='wireframe',opacity=1.0,color=(1.0,0.0,0.0))
        if show_plot == True and with_mayavi == False:
            autoscale_equal_axis(ax)
            # ax.autoscale()
            # ax.set_xlim(axmin,axmax)
            # ax.set_ylim(axmin,axmax)
            # ax.set_zlim(axmin,axmax)
        if show_plot==True:
            plt.show()
            return None 
        else:
            return ax

    
def plot_scatter_3d_plotly(m_sc_out):
    from plotly import __version__
    import plotly.graph_objs as go
    from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
    v_theta = m_sc_out[:,0,0]
    v_phi = m_sc_out[0,:,1]

    theta = m_sc_out[:,:,1].real
    phi = m_sc_out[:,:,0].real
    r = m_sc_out[:,:,6].real

    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.cos(phi)
    z = r * np.sin(phi) * np.sin(theta)


    surface = go.Surface(x=x, y=y, z=z, colorscale='Viridis',surfacecolor=r)
    data = [surface]
    layout = go.Layout(
        title='Dimer scattering',
        scene=dict(
            xaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            yaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            zaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            )
        )
    )
    fig = go.Figure(data=data, layout=layout)
    plot(fig, filename='parametric-plot-viridis') 

def plot_scatter_3d_color_sphere(m_sc_out,ax=None,color=None,show_plot=False):
    from matplotlib import cm
    v_theta = m_sc_out[:,0,0]
    v_phi = m_sc_out[0,:,1]

    theta = m_sc_out[:,:,1].real
    phi = m_sc_out[:,:,0].real
    r = np.log(m_sc_out[:,:,6].real)
    
    x = 1 * np.sin(phi) * np.cos(theta)
    y = 1 * np.cos(phi)
    z = 1 * np.sin(phi) * np.sin(theta)
    c = np.exp(r * np.sqrt(x**2+y**2+z**2))
    c = c - np.min(c)
    c = c/np.max(c)

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
    if color is None:

        ax.plot_surface(x,y,z,rstride=5, cstride=5, facecolors=cm.viridis(c),)
    else:
        ax.plot_surface(x,y,z,rstride=5, cstride=5, facecolors=cm.viridis(c), alpha=0.9, linewidth=1)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    autoscale_equal_axis(ax)
    if show_plot==True:
        plt.show()
        return None 
    else:
        return ax

def plot_scatter_3d_wireframe(m_sc_out,ax=None,show_plot=False,color=None,label="None",with_mayavi=False):

    
    v_theta = m_sc_out[:,0,0]
    v_phi = m_sc_out[0,:,1]

    theta = m_sc_out[:,:,1].real
    phi = m_sc_out[:,:,0].real
    r = m_sc_out[:,:,6].real
    
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.cos(phi)
    z = r * np.sin(phi) * np.sin(theta)


    if with_mayavi == False:

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            
        if color is None:
            ax.plot_wireframe(x,y,z,alpha=0.05)
        else:
            ax.plot_wireframe(x,y,z,alpha=0.05,color=color,label=label)
        autoscale_equal_axis(ax)
        if show_plot==True:
            plt.show()
            return None 
        else:
            return ax

    elif with_mayavi == True:
        mlab.gcf()
        mlab.mesh(x, y, z, representation='wireframe',opacity=0.05,color=tuple(color))
        # mlab.show()
    # fig = go.Figure(data=data, layout=layout)
    # plot(fig, filename='parametric-plot-viridis') 

if __name__ == "__main__":
    N = 2
    positions = np.asarray([[-10,0,0],[10,0,0]])
    # r = np.random.rand(N)/10.0
    r = [20,10]
    print r
    plot_3d_positions(positions,r)
    # mlab.show()


# [phi,theta] = np.mgrid[0:2*np.pi:12j,0:np.pi:12j]
# x = np.cos(phi)*np.sin(theta)
# y = np.sin(phi)*np.sin(theta)
# z = np.cos(theta)

# def plot_sphere(p):
#     r,a,b,c = p
#     return mlab.mesh(r*x+a, r*y+b, r*z+c)  

# for k in range(200):
#     c = np.random.rand(4)
#     c[0] /= 10.
#     plot_sphere(c)

# mlab.show()