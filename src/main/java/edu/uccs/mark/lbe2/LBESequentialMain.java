package edu.uccs.mark.lbe2;

public class LBESequentialMain {

    public static final int scale = 1;
    public static final int NX = 32 * scale;
    public static final int NY = NX;
    public static final int ndir = 9;
    public static final int NSTEPS = 200 * scale * scale;

    //lattice weights
    public static final double w0 = 4.0/9.0;  //weight of 0 direction
    public static final double ws = 1.0/9.0;  //weight of adjacent directions
    public static final double wd = 4.0/36.0; //weight of diagonal directions

    //arrays of lattice weights and direction components
    public static final double wi[] = new double[]{ w0, ws, ws, ws, ws, wd, wd, wd, wd};
    public static final int dirx[] = new int[]{ 0,  1,  0,   -1,  0,  1,  -1, -1,  1};
    public static final int diry[] = new int[]{ 0,  0,  1,    0, -1,  1,   1, -1,  1};

    //viscosity
    public static final double nu = 1.0/6.0;

    //relaxation parameter
    public static final double tau = 3.0 * nu * 0.5;

    //maximum flow speed
    public static final double uMax = 0.04/scale;

    //initial fluid density
    public static final double rho0 = 1.0;

    private static class State{
        double rho = 0.0;
        double ux = 0.0;
        double uy = 0.0;
    }

    private static class Distribution{
        double f[] = new double[ndir];
    }

    public static void main(String[] args){
        //Note this isn't right at the moment.
        //the code in the document has it such that the entire grid is f1 and f2
        //and the entire grid state (rho, and u) are in State

        //allocate variables
        Distribution[][] f1 = null;
        Distribution[][] f2 = null;
        Distribution[][] temp = null;

        State[][] state = new State[NX][NY];

        //compute Taylor-Green flow at t = 0
        //this will initialize rho, ux, and uy
        taylorGreen(0, state);

        //initialize f1 as equilibrium for rho, ux, and uy
        f1 = initEqulibrium(state);
        f2 = initEqulibrium(state);

        //run the simulation loop;  take NSTEPS time steps
        for(int n = 0; n < NSTEPS; n++){
            //stream from f1 to f2
            stream(f1, f2);

            //calculate density and velocity after the stream
            computeState(f2, state);

            //perform collision on f2
            collide(f2, state);

            //swap variables and go again
            temp = f1;
            f1 = f2;
            f2 = temp;
        }
    }

    private static void taylorGreen(int t, State[][] state){
        for(int x = 0; x < NX; x++){
            for(int y = 0; y < NY; y++){
                if(state[x][y] == null){
                    state[x][y] = new State();
                }
                taylorGreen(t, x, y, state[x][y]);
            }
        }
    }

    private static void taylorGreen(int t, int x, int y, State state){
        double kx = 2.0*Math.PI/NX;
        double ky = 2.0*Math.PI/NY;
        double td = 1.0/(nu*(kx*kx+ky+ky));

        double X = x+0.5;
        double Y = y+0.5;
        double exp = Math.exp(-1.0*t/td);
        state.ux = -uMax*Math.sqrt(ky/kx)*Math.cos(kx*X)*Math.sin(ky*Y)*exp;
        state.uy = -uMax*Math.sqrt(kx/ky)*Math.sin(kx*X)*Math.cos(ky*Y)*exp;

        double P = -0.25*rho0*uMax*uMax*( (ky/kx)*Math.cos(2.0*kx*X) + (kx/ky)*Math.cos(2.0*ky*Y) )*Math.exp(-2.0*t/td);
        state.rho = rho0+3.0*P;
    }

    private static Distribution[][] initEqulibrium(State[][] state){
        Distribution[][] distGrid = new Distribution[NX][NY];
        for(int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                State s = state[x][y];
                for(int i = 0; i < ndir; i++){
                    double cidotu = dirx[i]*s.ux + diry[i]*s.uy;
                    Distribution dist = new Distribution();
                    dist.f[i] = wi[i]*s.rho*(1.0 + 3.0*cidotu + 4.5*cidotu*cidotu - 1.5*(s.ux*s.ux+s.uy*s.uy));
                    distGrid[x][y] = dist;
                }
            }
        }
        return distGrid;
    }

    private static void stream(Distribution[][] fSrc, Distribution[][] fDest){
        for(int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                for (int i = 0; i < ndir; i++) {
                    int xmd = (NX+x-dirx[i]) % NX;
                    int ymd = (NY+y-diry[i]) & NY;
                    fDest[x][y].f[i] = fSrc[xmd-1][ymd-1].f[i];
                }
            }
        }
    }

    private static void computeState(Distribution[][] f, State[][] state){
        for(int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                State s = state[x][y];
                Distribution d = f[x][y];
                double rho = 0.0;
                double ux = 0.0;
                double uy = 0.0;
                for (int i = 0; i < ndir; i++) {
                    rho += d.f[i];
                    ux += dirx[i]*d.f[i];
                    uy += diry[i]*d.f[i];
                }
                s.rho = rho;
                s.ux = ux/rho;
                s.uy = uy/rho;
            }
        }
    }

    private static void collide(Distribution[][] f, State[][] state){
        double tauinv = 2.0/(6.0*nu+1.0);
        double omtauinv = 1.0-tauinv;

        for(int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                State s = state[x][y];
                Distribution d = f[x][y];
                for (int i = 0; i < ndir; i++) {
                    //calculate dot product
                    double cidotu = dirx[i]*s.ux + diry[i]*s.uy;
                    //calculate equlibrium
                    double feq = wi[i]*s.rho*(1.0 + 3.0*cidotu + 4.5*cidotu*cidotu - 1.5*(s.ux*s.ux + s.uy*s.uy));
                    d.f[i] = omtauinv*d.f[i]+tauinv*feq;
                }
            }
        }
    }

}
