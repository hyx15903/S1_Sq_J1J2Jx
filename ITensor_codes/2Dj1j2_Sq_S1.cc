#include "itensor/all.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#define M_PI 3.14159265358979323846 // pi

using namespace itensor;
using namespace std;

void make_ampo(auto ampo, auto lattice_type, auto j1, auto j2, auto jk, auto jz, auto Nx, auto Ny, auto N, auto yperiodic, auto args)
{
    if (lattice_type == "Square")
    {
        auto lattice_nn = squareNextNeighbor(Nx, Ny, args);
        for (auto &bnd : lattice_nn)
        {
            if (bnd.type == "1")
            {
                auto s1 = bnd.s1,
                     s2 = bnd.s2;
                int x1 = (s1 - 1) / Ny + 1;
                int x2 = (s2 - 1) / Ny + 1;
                ampo += (0.5 * j1), "S+", s1, "S-", s2;
                ampo += (0.5 * j1), "S-", s1, "S+", s2;
                ampo += (j1 * jz), "Sz", s1, "Sz", s2;
            }
            else if (bnd.type == "2")
            {
                auto s1 = bnd.s1,
                     s2 = bnd.s2;
                int x1 = (s1 - 1) / Ny + 1;
                int x2 = (s2 - 1) / Ny + 1;
                ampo += (0.5 * j2), "S+", s1, "S-", s2;
                ampo += (0.5 * j2), "S-", s1, "S+", s2;
                ampo += (j2 * jz), "Sz", s1, "Sz", s2;
            }
        }
        if (jk != 0)
        {
            // all clockwise
            //-i/2 S_z S_- S_+
            // i/2 S_z S_+ S_-
            //-i/2 S_+ S_z S_-
            // i/2 S_- S_z S_+
            //-i/2 S_- S_+ S_z
            // i/2 S_+ S_- S_z
            jk = jk / 2.;
            for (int n = 1; n <= N; n++)
            {
                int x = (n - 1) / Ny + 1;
                int y = (n - 1) % Ny + 1;
                // x+45 direction
                if ((y < Ny) && (x < Nx))
                {
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S-", n + 1, "S+", n + Ny + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S+", n + 1, "S-", n + Ny + 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S+", n, "Sz", n + 1, "S-", n + Ny + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S-", n, "Sz", n + 1, "S+", n + Ny + 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S-", n, "S+", n + 1, "Sz", n + Ny + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S+", n, "S-", n + 1, "Sz", n + Ny + 1;
                }
                if ((yperiodic) && (y == Ny) && (x < Nx))
                {
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S-", n - Ny + 1, "S+", n + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S+", n - Ny + 1, "S-", n + 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S+", n, "Sz", n - Ny + 1, "S-", n + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S-", n, "Sz", n - Ny + 1, "S+", n + 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S-", n, "S+", n - Ny + 1, "Sz", n + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S+", n, "S-", n - Ny + 1, "Sz", n + 1;
                }
                // x-45 direction
                if ((y > 1) && (x < Nx))
                {
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S-", n + Ny, "S+", n + Ny - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S+", n + Ny, "S-", n + Ny - 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S+", n, "Sz", n + Ny, "S-", n + Ny - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S-", n, "Sz", n + Ny, "S+", n + Ny - 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S-", n, "S+", n + Ny, "Sz", n + Ny - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S+", n, "S-", n + Ny, "Sz", n + Ny - 1;
                }
                if ((yperiodic) && (y == 1) && (x < Nx))
                {
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S-", n + Ny, "S+", n + 2 * Ny - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S+", n + Ny, "S-", n + 2 * Ny - 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S+", n, "Sz", n + Ny, "S-", n + 2 * Ny - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S-", n, "Sz", n + Ny, "S+", n + 2 * Ny - 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S-", n, "S+", n + Ny, "Sz", n + 2 * Ny - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S+", n, "S-", n + Ny, "Sz", n + 2 * Ny - 1;
                }
                // x+135 direction
                if ((y < Ny) && (x > 1))
                {
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S-", n - Ny, "S+", n - Ny + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S+", n - Ny, "S-", n - Ny + 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S+", n, "Sz", n - Ny, "S-", n - Ny + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S-", n, "Sz", n - Ny, "S+", n - Ny + 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S-", n, "S+", n - Ny, "Sz", n - Ny + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S+", n, "S-", n - Ny, "Sz", n - Ny + 1;
                }
                if ((yperiodic) && (y == Ny) && (x > 1))
                {
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S-", n - Ny, "S+", n - 2 * Ny + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S+", n - Ny, "S-", n - 2 * Ny + 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S+", n, "Sz", n - Ny, "S-", n - 2 * Ny + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S-", n, "Sz", n - Ny, "S+", n - 2 * Ny + 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S-", n, "S+", n - Ny, "Sz", n - 2 * Ny + 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S+", n, "S-", n - Ny, "Sz", n - 2 * Ny + 1;
                }
                // x-135 direction
                if ((y > 1) && (x > 1))
                {
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S-", n - 1, "S+", n - Ny - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S+", n - 1, "S-", n - Ny - 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S+", n, "Sz", n - 1, "S-", n - Ny - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S-", n, "Sz", n - 1, "S+", n - Ny - 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S-", n, "S+", n - 1, "Sz", n - Ny - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S+", n, "S-", n - 1, "Sz", n - Ny - 1;
                }
                if ((yperiodic) && (y == 1) && (x > 1))
                {
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S-", n + Ny - 1, "S+", n - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S+", n + Ny - 1, "S-", n - 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S+", n, "Sz", n + Ny - 1, "S-", n - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S-", n, "Sz", n + Ny - 1, "S+", n - 1;
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S-", n, "S+", n + Ny - 1, "Sz", n - 1;
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S+", n, "S-", n + Ny - 1, "Sz", n - 1;
                }
            }
            jk *= 2.;
        }
    }
    else
        Error("Wrong lattice type");
}

int main(int argc, char *argv[])
{
    // H = (jz * j1 S_z^i * S_z^i+1 + j1 /2 S_+^i * S_-^i+1 + j1 /2 S_-^i * S_+^i+1)
    //   + (jz * j2 S_z^i * S_z^i+2 + j2 /2 S_+^i * S_-^i+2 + j2 /2 S_-^i * S_+^i+2)
    // Sweeps
    // Parse the input file
    //
    if (argc < 2)
    {
        printfln("Usage: %s SweepTable2Dj1j2", argv[0]);
        return 0;
    }
    auto input = InputGroup(argv[1], "input");

    // Read in individual parameters from the input file
    //
    //
    // Other parameters
    auto Nx = input.getInt("Nx");
    auto Ny = input.getInt("Ny");
    auto j1 = input.getReal("j1");      // nn exchange
    auto j2 = input.getReal("j2");      // nnn exchange, for kagome it is nnn and nnnn
    auto jk = input.getReal("jk");      // chiral interaction
    auto jz = input.getReal("jz");      // unisotropy in z direction
    auto totSz = input.getInt("TotSz"); // ground state total Sz
    int N = Nx * Ny;
    auto nsweeps = input.getInt("nsweeps");
    auto quiet = input.getYesNo("quiet", true);
    auto yperiodic = input.getYesNo("periodic", false); // periodic boundary in y direction
    auto lattice_type = input.getString("lattice_type", "Square");
    auto ErrorGoal = input.getReal("ErrorGoal");
    auto DoDMRG = input.getYesNo("DoDMRG", false);
    auto SavePsiSites = input.getYesNo("SavePsiSites", true);
    auto DoSpinStructure = input.getYesNo("DoSpinStructure", false);     // spin structure
    auto DoCorrelation = input.getYesNo("DoCorrelation", false);         // spin correlations
    auto TypeofCorrelation = input.getString("TypeofCorrelation", "S");  // S, Sz, Sxy
    auto CmtIn = input.getString("CommentsInPsiSites", "testin");        // comment for input file
    auto CmtOut = input.getString("CommentsOutPsiSitesCorr", "testout"); // comment for output file

    auto table = InputGroup(input, "sweeps");
    // Create the sweeps class & print
    //
    auto sweeps = Sweeps(nsweeps, table);
    println(sweeps);
    // set variables to string type
    //
    std::ostringstream strNxs;
    strNxs << Nx;
    string strNx = strNxs.str();
    std::ostringstream strNys;
    strNys << Ny;
    string strNy = strNys.str();
    std::ostringstream strj2s;
    strj2s << j2;
    std::string strj2 = strj2s.str();
    std::ostringstream strjzs;
    strjzs << jz;
    std::string strjz = strjzs.str();
    std::ostringstream strjks;
    strjks << jk;
    std::string strjk = strjks.str();
    std::ostringstream strtotSzs;
    strtotSzs << totSz;
    std::string strtotSz = strtotSzs.str();

    Args args;
    args.add("Quiet", quiet);
    args.add("YPeriodic", yperiodic);
    args.add("maxdim", 500);
    args.add("Cutoff", 1E-8);
    println(args);
    // Define space
    //
    SiteSet sites;
    sites = SpinOne(N);
    // Initialize the site degrees of freedom.
    //
    cout << "MPO" << endl;
    auto ampo = AutoMPO(sites);
    make_ampo(ampo, lattice_type, j1, j2, jk, jz, Nx, Ny, N, yperiodic, args);
    auto H = toMPO(ampo);
    cout << "MPS" << endl;
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    MPS psi(sites);
    auto state = InitState(sites);
    for (int i = 1; i <= N; ++i)
    {
        if (i % 2 == 1)
            state.set(i, "Dn");
        else
            state.set(i, "Up");
    }
    // totSz is not zero
    for (int i = 1; i <= totSz; ++i)
    {
        state.set(N / 2 - i, "Up");
        state.set(N / 2 + i - 1, "Up");
    }

    psi = MPS(state);
    auto psi1 = psi;
    //
    // complex initial state?
    psi *= 2.0 * Cplx_i;
    psi = sum(psi, psi1, args);
    psi.normalize();
    println(totalQN(psi));
    printfln("Initial energy = %.8f", innerC(psi, H, psi));

    if (DoDMRG)
    {
        auto obs = DMRGObserver(psi, {"EnergyErrgoal", ErrorGoal});
        auto [energy, psi0] = dmrg(H, psi, sweeps, obs);
        psi = psi0;
        // cout << "doing excited states" << endl;
        // auto wfs = std::vector<MPS>(1);
        // wfs.at(0) = psi;
        //
        //  Here the Weight option sets the energy penalty for
        //  psi1 having any overlap with psi
        //
        // auto [en1, psi1] = dmrg(H, wfs, randomMPS(state), sweeps, {"Quiet=", true, "Weight=", 20.0});
        // psi = psi1;
    }

    // save the ground state to disk
    if (SavePsiSites)
    {
        writeToFile("sites_S1_" + lattice_type + "_Nx" + strNx + "Ny" + strNy + "j2" + strj2 + "jz" + strjz + "jk" + strjk + CmtOut, sites);
        writeToFile("psi_S1_" + lattice_type + "_Nx" + strNx + "Ny" + strNy + "j2" + strj2 + "jz" + strjz + "jk" + strjk + CmtOut, psi);
    }

    // Print the final energy reported by DMRG
    //
    // printfln("\nGround State Energy = %.10f", energy);
    printfln("\nUsing inner = %.10f", innerC(psi, H, psi));

    println("\nTotal QN of Ground State = ", totalQN(psi));
    cout << endl;

    int x1, x2, posi;
    double sumCorr;
    posi = 1;
    psi.position(posi);
    auto C = sites.op("Sz", posi) * psi(posi);
    if (DoCorrelation)
    {
        // Static correlation function
        cout << "Definition: Correlation= <S_i*S_j>" << endl;
        int CorrRange = 8 * Ny + 2;      // range
        int CorrStart = 4 * Ny + Ny / 2; // first site
        int CorrEnd = 4 * Ny + Ny / 2;   // last site
        double CorrAve[CorrRange];
        for (int i = 0; i <= (CorrRange - 1); i++)
        {
            CorrAve[i] = 0.0;
        }
        ofstream File;
        File.open("Stat_" + TypeofCorrelation + "_Corr_S1_Nx" + strNx + "Ny" + strNy + "j2" + strj2 + "jz" + strjz + "jk" + strjk + CmtOut + ".txt");
        for (int j = CorrStart; j <= CorrEnd; j++)
        {
            cout << "Correlation at site: " << j << endl;
            for (int i = 1; i <= (CorrRange - 1); i++)
            {
                x1 = j;
                x2 = x1 + i;
                Sab = 0.0;
                sumCorr = 0.0;
                // single site magnetization
                /*if ((TypeofCorrelation == "S") || (TypeofCorrelation == "Sz"))
                {
                    posi = x1;
                    psi.position(posi);
                    C = dag(prime(psi(posi), "Site")) * sites.op("Sz", posi) * psi(posi);
                    cout << "<Sz(" << posi << ")> = " << C.cplx().real() << endl;
                    Sab += C.cplx().real();

                    posi = x2;
                    psi.position(posi);
                    C = dag(prime(psi(posi), "Site")) * sites.op("Sz", posi) * psi(posi);
                    //printfln("\nSz = %.10f", C);
                    cout << "<Sz(" << posi << ")> = " << C.cplx().real() << endl;
                    Sab = Sab * C.cplx().real();
                }*/

                // Sz * Sz
                if ((TypeofCorrelation == "S") || (TypeofCorrelation == "Sz"))
                {
                    auto op_az = sites.op("Sz", x1);
                    auto op_bz = sites.op("Sz", x2);
                    psi.position(x1);
                    auto ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * op_az * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < x2; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(x2);
                    C *= op_bz;
                    auto jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    // cout << "<Sz(1) * Sz(L/2)> = " << (eltC(C).real()) << endl;
                    // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                    //
                    sumCorr += C.cplx().real();
                }

                // S+ * S-
                if ((TypeofCorrelation == "S") || (TypeofCorrelation == "Sxy"))
                {
                    auto op_ap = sites.op("S+", x1);
                    auto op_bm = sites.op("S-", x2);
                    psi.position(x1);
                    auto ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * op_ap * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < x2; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(x2);
                    C *= op_bm;
                    auto jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    // cout << "0.5<S+(1) * S-(L/2)> = " << (eltC(C).real() / 2.0) << endl;
                    // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                    //
                    sumCorr += C.cplx().real() / 2.0;

                    // S- * S+
                    auto op_am = sites.op("S-", x1);
                    auto op_bp = sites.op("S+", x2);
                    psi.position(x1);
                    ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * op_am * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < x2; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(x2);
                    C *= op_bp;
                    jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    // cout << "0.5<S-(1) * S+(L/2)> = " << ( eltC(C).real() / 2.0) << endl;
                    // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                    //
                    sumCorr += C.cplx().real() / 2.0;
                }
                // sumCorr -= Sab;
                cout << "Correlation of [" << TypeofCorrelation << "] (" << x1 << ")(" << x2 << ") = " << sumCorr << endl;
                File << i << " " << sumCorr << endl;
                CorrAve[i] += sumCorr;
            }
        }
        cout << "Averagr of Correlation of spin: " << endl;
        for (int i = 1; i <= (CorrRange - 1); i++)
        {
            cout << i << " " << CorrAve[i] / (CorrEnd - CorrStart + 1) << endl;
        }
        File.close();
    }

    if (DoSpinStructure)
    {
        // correlations within Ny * Ny lattice
        cout << "correlations within Ny * Ny lattice" << endl;
        ofstream SpinStru1;
        ofstream SpinStruDataReal1;
        SpinStru1.open("SpinStructure_Ny_X_Ny_S1_Nx" + strNx + "Ny" + strNy + "j2" + strj2 + "jz" + strjz + "jk" + strjk + CmtOut + ".txt");
        SpinStruDataReal1.open("SpinStru_RealSpace_Ny_X_Ny_S1_Nx" + strNx + "Ny" + strNy + "j2" + strj2 + "jz" + strjz + "jk" + strjk + CmtOut + ".txt");
        vector<vector<double>> SpinCorrValue1(Ny * Ny);
        for (int i = 0; i < (Ny * Ny); i++)
        {
            cout << "Calculating point (" << i << ")" << endl;
            for (int j = 0; j < (Ny * Ny); j++)
            {
                x1 = i + 1 + Ny * (Nx - Ny) / 2;
                x2 = j + 1 + Ny * (Nx - Ny) / 2;
                if (x2 < x1)
                {
                    int x3 = x2;
                    x2 = x1;
                    x1 = x3;
                }
                else if (x2 == x1)
                {
                    x1 = x1 - 1;
                }
                // Sz * Sz
                auto op_az = sites.op("Sz", x1);
                auto op_bz = sites.op("Sz", x2);
                psi.position(x1);
                auto ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                C = psi(x1) * op_az * dag(prime(prime(psi(x1), "Site"), ir));
                for (int k = x1 + 1; k < x2; ++k)
                {
                    C *= psi(k);
                    C *= dag(prime(psi(k), "Link"));
                }
                C *= psi(x2);
                C *= op_bz;
                auto jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                C *= dag(prime(prime(psi(x2), jl), "Site"));
                // cout << "<Sz(1) * Sz(L/2)> = " << (C.cplx().real()) << endl;
                // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                //
                sumCorr = C.cplx().real();
                // S+ * S-
                auto op_ap = sites.op("S+", x1);
                auto op_bm = sites.op("S-", x2);
                psi.position(x1);
                ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                C = psi(x1) * op_ap * dag(prime(prime(psi(x1), "Site"), ir));
                for (int k = x1 + 1; k < x2; ++k)
                {
                    C *= psi(k);
                    C *= dag(prime(psi(k), "Link"));
                }
                C *= psi(x2);
                C *= op_bm;
                jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                C *= dag(prime(prime(psi(x2), jl), "Site"));
                // cout << "0.5<S+(1) * S-(L/2)> = " << (C.cplx().real() / 2.0) << endl;
                // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                //
                sumCorr += C.cplx().real() / 2.0;
                // S- * S+
                auto op_am = sites.op("S-", x1);
                auto op_bp = sites.op("S+", x2);
                psi.position(x1);
                ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                C = psi(x1) * op_am * dag(prime(prime(psi(x1), "Site"), ir));
                for (int k = x1 + 1; k < x2; ++k)
                {
                    C *= psi(k);
                    C *= dag(prime(psi(k), "Link"));
                }
                C *= psi(x2);
                C *= op_bp;
                jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                C *= dag(prime(prime(psi(x2), jl), "Site"));
                // cout << "0.5<S-(1) * S+(L/2)> = " << ( C.cplx().real() / 2.0) << endl;
                // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                //
                sumCorr += C.cplx().real() / 2.0;
                // New definition: Correlation= <S_i*S_j>
                // sumCorr -= Sab;
                if (i == j)
                {
                    SpinCorrValue1[i].emplace_back(2.);
                    SpinStruDataReal1 << i << " " << j << " " << 2. << endl;
                }
                else
                {
                    SpinCorrValue1[i].emplace_back(sumCorr);
                    SpinStruDataReal1 << i << " " << j << " " << sumCorr << endl;
                }
            }
            SpinStruDataReal1 << endl;
        }
        SpinStruDataReal1.close();

        double Sq = 0.0;
        for (double qx = (-M_PI * 2.); qx < (M_PI * 2. + 0.02); qx += (2. * M_PI / Ny))
        {
            for (double qy = ((-2.) * M_PI); qy < (2. * M_PI + 0.02); qy += (2. * M_PI / Ny))
            {
                // spin structure within Ny * Ny lattice
                double Sq1 = 0.0;
                for (int i = 0; i < (Ny * Ny); i++)
                {
                    for (int j = 0; j < (int)SpinCorrValue1[i].size(); j++)
                    {
                        int xi = i / Ny;
                        int yi = i % Ny;
                        int xj = j / Ny;
                        int yj = j % Ny;
                        Sq1 += SpinCorrValue1[i][j] * cos(qx * (xj - xi) + qy * (yj - yi));
                    }
                }
                Sq1 /= (Ny * Ny);
                SpinStru1 << qx << " " << qy << " " << Sq1 << endl;
            }
        }
        SpinStru1.close();
    }

    // print parameters
    //
    cout << " j1 = " << j1 << endl;
    cout << " j2 = " << j2 << endl;
    cout << " jz = " << jz << endl;
    cout << " Nx = " << Nx << endl;
    cout << " Ny = " << Ny << endl;
    cout << " jk = " << jk << endl;
    cout << " totSz = " << totSz << endl;
    cout << " lattice type = " << lattice_type << endl;
    // double E0 = inner(psi, H, psi);
    // double E02 = inner(psi, H, H, psi);
    // cout << " Ground State variance = " << E02 - E0 * E0 << endl;
    return 0;
}
