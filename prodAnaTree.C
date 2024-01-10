                                                                                  
#include <TSystem>
class StMaker;
class StChain;
class StPicoDstMaker;

StChain *chain;
void prodAnaTree(const Char_t *inputFile="test.list", const Char_t *outputFile="test.root")
{
        //Load all the System libraries
        gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
        loadSharedLibraries();

        gSystem->Load("StChain");
        gSystem->Load("StPicoEvent");
        gSystem->Load("StPicoDstMaker");
        gSystem->Load("StPico2AnaMaker");

        //    gSystem->Load("StBichsel");
        //    gSystem->Load("StDbLib");
        //    gSystem->Load("StDbBroker");
        //    gSystem->Load("StDetectorDbMaker");
        //    gSystem->Load("St_db_Maker");
        //    gSystem->Load("StDetectorDbMaker.so");
        //    gSystem->Load("StTpcDb");
        //    gSystem->Load("StMagF");
        //    gSystem->Load("StDaqLib");
        //    gSystem->Load("libgen_Tables");
        //    gSystem->Load("libsim_Tables");
        //    gSystem->Load("libglobal_Tables");

        chain = new StChain();

        StPicoDstMaker *picoMaker = new StPicoDstMaker(2, inputFile, "picoDst");

        //St_db_Maker *dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");
        //dbMk->SetDateTime(2016*1e4+101,0);
        //StMagFMaker *magfMk = new StMagFMaker;

        StPico2AnaMaker *anaMaker = new StPico2AnaMaker("ana",picoMaker,outputFile);

        //-----------------------------------------
        anaMaker->addTrigger(350007, 0); //UPC_main
        anaMaker->addTrigger(350017, 0); //UPC_main

        anaMaker->addTrigger(350007, 1); //UPC_main1
        anaMaker->addTrigger(350017, 2); //UPC_main2

        chain->Init();
        cout << "chain->Init();" << endl;

        int nEvents = picoMaker->chain()->GetEntries();
        cout << " Total entries = " << nEvents << endl;

        int total = 0;
        for ( int i = 0; i < nEvents; i++ )
        {
                if ( i%1000 == 0 )
                {
                        cout << "Working on eventNumber " << i << " of "<< nEvents << " (" << fixed <<setprecision(1) << "finished "<< ((float)i/(float)nEvents)*100.$
                }

                chain->Clear();
                int iret = chain->Make(i);

                if ( iret )
                {
                        cout << "Bad return code!" << iret << endl;
                        break;
                }

                total++;
        }

	cout << "****************************************** " << endl;
        cout << "Work done... now its time to close up shop!" << endl;
        cout << "****************************************** " << endl;
        chain->Finish();
        cout << "****************************************** " << endl;
        cout << "total number of events  " << total << endl;
        cout << "****************************************** " << endl;

        delete chain;
}
