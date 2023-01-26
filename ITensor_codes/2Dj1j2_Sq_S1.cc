#include "itensor/all.h"
#include "spink.h"
#include "spin_site_chiral.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#define M_PI 3.14159265358979323846 // pi

using namespace itensor;
using namespace std;

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
    auto j1 = input.getReal("j1");       // nn exchange
    auto j2 = input.getReal("j2");       // nnn exchange, for kagome it is nnn and nnnn
    auto jk = input.getReal("jk");       // chiral interaction
    auto totSz = input.getInt("TotSz");
    auto Bz = input.getReal("Bz");               // Zeeman field for the triangularXC_Ky
    auto hz = input.getReal("hz");               // edge magnetic pining field
    auto torus = input.getYesNo("torus", false); // torus
    int N = Nx * Ny;
    auto nsweeps = input.getInt("nsweeps");
    auto quiet = input.getYesNo("quiet", true);
    auto yperiodic = input.getYesNo("periodic", false);
    auto lattice_type = input.getString("lattice_type", "triangularXC");
    auto ErrorGoal = input.getReal("ErrorGoal");
    auto DoDMRG = input.getYesNo("DoDMRG", false);
    auto DoSpinPumping = input.getYesNo("DoSpinPumping", false);
    auto SavePsiSites = input.getYesNo("SavePsiSites", true);
    auto DoSpinStructure = input.getYesNo("DoSpinStructure", false);
    auto DoCorrelation = input.getYesNo("DoCorrelation", false);
    auto TypeofCorrelation = input.getString("TypeofCorrelation", "S"); // S, Sz, Sxy
    auto CmtIn = input.getString("CommentsInPsiSites", "testin");
    auto CmtOut = input.getString("CommentsOutPsiSitesCorr", "testout");

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
    std::ostringstream stromegas;
    stromegas << omega;
    std::string stromega = stromegas.str();
    std::ostringstream strtotKys;
    strtotKys << totKy;
    std::string strtotKy = strtotKys.str();
    std::ostringstream strtotSzs;
    strtotSzs << totSz;
    std::string strtotSz = strtotSzs.str();

    Args args;
    args.add("Quiet", quiet);
    args.add("YPeriodic", yperiodic);
    args.add("maxdim", 500);
    args.add("Cutoff", 1E-8);
    if (lattice_type == "triangularXC_Ky")
    {
        args.add("ConserveQNs", true);
        args.add("ConserveK", true);
        args.add("ConserveSz", true);
        args.add("Kmod", Ny);
    }
    println(args);

    // Define space
    //
    SiteSet sites;
    sites = SpinOne(N);

    // Initialize the site degrees of freedom.
    //
    cout << "MPO" << endl;
    auto ampo = AutoMPO(sites);
    auto ampo_mid = AutoMPO(sites);
    if (!DoSpinPumping)
    {
        make_ampo(ampo, ampo_mid, lattice_type, sites, j1, j2, jk, jz, omega, jkedge, hz, Bz, torus, Nx, Ny, strNx, strNy, N, yperiodic, args, DoMidEnergyPerSite, DoMidDMRGSweeps, mid_from, mid_to);
    }
    auto H = toMPO(ampo);
    cout << "MPS" << endl;
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    MPS psi;
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
    
    psi.position(N / 2);
    psi.normalize();
    println(totalQN(psi));
    printfln("Initial energy = %.8f", innerC(psi, H, psi));

    if (DoDMRG)
    {
        cout << "dmrg" << endl;
        auto obs = DMRGObserver(psi, {"EnergyErrgoal", ErrorGoal});
        auto [energy, psi0] = dmrg(H, psi, sweeps, obs, args1);
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
        cout << "dmrg done" << endl;
    }
    //
    // spin pumping
    //
    if (DoSpinPumping)
    {
        for (int iter = 1; iter <= 8; iter++)
        {
            auto ampo_pump = AutoMPO(sites);
            auto ampo_pump_mid = AutoMPO(sites);
            make_ampo(ampo_pump, ampo_pump_mid, lattice_type, sites, j1, j2, jk, jz, iter / 8., jkedge, hz, Bz, torus, Nx, Ny, strNx, strNy, N, yperiodic, args, DoMidEnergyPerSite, DoMidDMRGSweeps, mid_from, mid_to);
            auto H_pump = toMPO(ampo_pump);
            auto obs = DMRGObserver(psi, {"EnergyErrgoal", ErrorGoal});
            auto [energy, psi0] = dmrg(H_pump, psi, sweeps, obs, args1);
            save_psisite(Spin, lattice_type, strNx, strNy, strj2, strjz, strjk, strtotKy, stromega, strtotSz, CmtOut, sites, psi0);
            psi = psi0;
        }
    }

    // save the ground state to disk
    if (SavePsiSites)
    {
        save_psisite(Spin, lattice_type, strNx, strNy, strj2, strjz, strjk, strtotKy, stromega, strtotSz, CmtOut, sites, psi);
    }

    // Print the final energy reported by DMRG
    //
    // printfln("\nGround State Energy = %.10f", energy);
    printfln("\nUsing inner = %.10f", innerC(psi, H, psi));
    
    println("\nTotal QN of Ground State = ", totalQN(psi));
    cout << endl;

    int x1, x2, posi;
    double sumCorr, Sab;
    posi = 1;
    psi.position(posi);
    auto C = sites.op("Sz", posi) * psi(posi);
    auto D = sites.op("Sz", posi) * psi(posi);
    auto EE = sites.op("Sz", posi) * psi(posi);
    auto FF = sites.op("Sz", posi) * psi(posi + 1);
    auto GG = sites.op("Sz", posi) * psi(posi);
    auto HH = sites.op("Sz", posi) * psi(posi);
    auto II = sites.op("Sz", posi) * psi(posi);
    auto JJ = sites.op("Sz", posi) * psi(posi);

    if (DoCorrelation)
    {
        // Static correlation function (<S(1,1) * S(L/2,L/2)> - <S(1,1)> * <S(L/2,L/2)>) of electron spin (old)
        cout << "New definition: Correlation= <S_i*S_j>" << endl;
        int CorrRange = 8 * Ny + 2;      // range
        int CorrStart = 4 * Ny + Ny / 2; // first site
        int CorrEnd = 4 * Ny + Ny / 2;   // last site
        double CorrAve[CorrRange];
        for (int i = 0; i <= (CorrRange - 1); i++)
        {
            CorrAve[i] = 0.0;
        }
        ofstream File;
        File.open("Stat" + TypeofCorrelation + "Corr_S1_Nx" + strNx + "Ny" + strNy + "j2" + strj2 + "jz" + strjz + "jk" + strjk + "omega" + stromega + "_totKy_" + strtotKy + CmtOut + ".txt");
        for (int j = CorrStart; j <= CorrEnd; j++)
        {
            cout << "Correlation at site: " << j << endl;
            for (int i = 1; i <= (CorrRange - 1); i++)
            {
                // correlation function (<S(1,1) * S(L/2,L/2)> - <S(1,1)> * <S(L/2,L/2)>) of electron spin
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
        cout << endl;
        int howmuchpoint = Ny * 5;
        // int xx = 84 / Ny; old version
        cout << "new xrange" << endl;
        cout << "New definition: Correlation= <S_i*S_j>" << endl;
        int xx = Nx / 2;
        int CorrRange = Ny * xx; // range
        // int CorrStart = (N - (CorrRange + howmuchpoint)) / 2 + 1; // first site
        int CorrStart = 1;
        // int CorrEnd = CorrStart + howmuchpoint - 1;			  // last site
        vector<vector<double>> SpinCorrValue(howmuchpoint);
        cout << "SpinStructure begin: " << endl;
        ofstream SpinStru;
        ofstream SpinStruDataReal;

        /*SpinStru.open("SpinStructure_S1_Nx" + strNx + "Ny" + strNy + "j2" + strj2 + "jz" + strjz + "jk" + strjk + "omega" + stromega + CmtOut + ".txt");
        SpinStruDataReal.open("SpinStruDataRealSpace_S1_Nx" + strNx + "Ny" + strNy + "j2" + strj2 + "jz" + strjz + "jk" + strjk + "omega" + stromega + CmtOut + ".txt");
        cout << "total points = " << howmuchpoint << endl;
        cout << "Correlation Range = " << CorrRange << endl;
        cout << "Correlation Starting point = " << CorrStart << endl;
        cout << endl;
        for (int i = 0; i < howmuchpoint; i++)
        {
            cout << "Calculating point (" << i << ")" << endl;
            for (int j = 1; j <= CorrRange; j++)
            {

                //correlation function
                x1 = i + CorrStart;
                x2 = j + x1;
                //single site magnetization
                posi = x1;
                psi.position(posi);
                C = dag(prime(psi(posi), "Site")) * sites.op("Sz", posi) * psi(posi);
                Sab = C.cplx().real();
                posi = x2;
                psi.position(posi);
                C = dag(prime(psi(posi), "Site")) * sites.op("Sz", posi) * psi(posi);
                //printfln("\nSz = %.10f", C);
                Sab = Sab * C.cplx().real();
                //Sz * Sz
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
                //cout << "<Sz(1) * Sz(L/2)> = " << (C.cplx().real()) << endl;
                //cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                //
                sumCorr = C.cplx().real();
                //S+ * S-
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
                //cout << "0.5<S+(1) * S-(L/2)> = " << (C.cplx().real() / 2.0) << endl;
                //cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                //
                sumCorr += C.cplx().real() / 2.0;
                //S- * S+
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
                //cout << "0.5<S-(1) * S+(L/2)> = " << ( C.cplx().real() / 2.0) << endl;
                //cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                //
                sumCorr += C.cplx().real() / 2.0;
                // New definition: Correlation= <S_i*S_j>
                //sumCorr -= Sab;
                SpinStruDataReal << j << " " << sumCorr << endl;
                SpinCorrValue[i].emplace_back(sumCorr);
            }
            SpinStruDataReal << endl;
            SpinStruDataReal << endl;
            SpinStruDataReal << endl;
            SpinStruDataReal << endl;
        }
        SpinStruDataReal.close();*/

        // correlations within Ny * Ny lattice
        cout << "correlations within Ny * Ny lattice" << endl;
        ofstream SpinStru1;
        ofstream SpinStruDataReal1;
        SpinStru1.open("SpinStructureNy_X_Ny_S1_Nx" + strNx + "Ny" + strNy + "j2" + strj2 + "jz" + strjz + "jk" + strjk + "omega" + stromega + CmtOut + ".txt");
        SpinStruDataReal1.open("SpinStruRealSpaceNy_X_Ny_S1_Nx" + strNx + "Ny" + strNy + "j2" + strj2 + "jz" + strjz + "jk" + strjk + "omega" + stromega + CmtOut + ".txt");
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
                    SpinCorrValue1[i].emplace_back(1.);
                    SpinStruDataReal1 << i << " " << j << " " << 1. << endl;
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
        if (lattice_type == "HoneycombXC")
        {
            for (double qx = (-M_PI * 2. / sqrt(3)); qx < (M_PI * 2. / sqrt(3) + 0.05); qx += (4. * M_PI / Ny / sqrt(3)))
            {
                for (double qy = (-4. * M_PI / 3.); qy < (4. * M_PI / 3. + 0.05); qy += (8. * M_PI / Ny / 3.))
                {
                    Sq = 0.0;
                    for (int i = 0; i < howmuchpoint; i++)
                    {
                        for (int j = 0; j < (int)SpinCorrValue[i].size(); j++)
                        {
                            int xx = (j + 1) / Ny;
                            int yy = (j + 1) % Ny;
                            if ((i + CorrStart) % 2 == 0)
                            {
                                if ((i + CorrStart + j + 1) % 2 == 0)
                                {
                                    Sq += 2. * SpinCorrValue[i][j] * cos(qx * (xx * sqrt(3) - (yy * sqrt(3) / 4)) + qy * yy / 2 * 1.5);
                                }
                                else
                                {
                                    Sq += 2. * SpinCorrValue[i][j] * cos(qx * (xx * sqrt(3) - ((yy - 1) * sqrt(3) / 4)) + qy * ((yy - 1) / 2 * 1.5 + 1.));
                                }
                            }
                            else
                            {
                                if ((i + CorrStart + j + 1) % 2 == 0)
                                {
                                    Sq += 2. * SpinCorrValue[i][j] * cos(qx * (xx * sqrt(3) - ((yy + 1) * sqrt(3) / 4)) + qy * ((yy - 1) / 2 * 1.5 + 0.5));
                                }
                                else
                                {
                                    Sq += 2. * SpinCorrValue[i][j] * cos(qx * (xx * sqrt(3) - (yy * sqrt(3) / 4)) + qy * yy / 2 * 1.5);
                                }
                            }
                        }
                    }
                    Sq /= howmuchpoint;
                    SpinStru << qx << " " << qy << " " << Sq << endl;
                }
            }
        }
        else if (lattice_type == "Square")
        {
            for (double qx = (-M_PI * 2.); qx < (M_PI * 2. + 0.02); qx += (2. * M_PI / Ny))
            {
                for (double qy = ((-2.) * M_PI); qy < (2. * M_PI + 0.02); qy += (2. * M_PI / Ny))
                {
                    /*Sq = 0.0;
                    for (int i = 0; i < howmuchpoint; i++)
                    {
                        for (int j = 0; j < (int)SpinCorrValue[i].size(); j++)
                        {
                            int xx = (j + 1) / Ny;
                            int yy = (j + 1) % Ny;
                            Sq += 2. * SpinCorrValue[i][j] * cos(qx * xx + qy * yy);
                        }
                    }
                    Sq /= howmuchpoint;
                    SpinStru << qx << " " << qy << " " << Sq << endl;*/

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
        }
        // SpinStru.close();
        SpinStru1.close();
    }

    // print parameters
    //
    cout << " j1 = " << j1 << endl;
    cout << " j2 = " << j2 << endl;
    cout << " jz = " << jz << endl;
    cout << " Nx = " << Nx << endl;
    cout << " Ny = " << Ny << endl;
    cout << " hz = " << hz << endl;
    cout << " jk = " << jk << endl;
    cout << " jkedge = " << jkedge << endl;
    cout << " totKy = " << totKy << endl;
    cout << " totSz = " << totSz << endl;
    cout << " lattice type = " << lattice_type << endl;

    // double E0 = inner(psi, H, psi);
    // double E02 = inner(psi, H, H, psi);
    // cout << " Ground State variance = " << E02 - E0 * E0 << endl;
    return 0;
}
