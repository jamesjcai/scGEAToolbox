function callback_RunEnrichr(src, ~, predefinedlist, enrichrtype, ...
    backgroundlist, outfiletag, wkdir)

if nargin < 7, wkdir = ''; end
if nargin < 6, outfiletag = ""; end
if nargin < 5
    askbackground = true;
    backgroundlist = []; 
else
    askbackground = false;
end
if nargin < 4, enrichrtype = []; end
if nargin < 3, predefinedlist = []; end

%    [FigureHandle, sce] = gui.gui_getfigsce(src);

    [FigureHandle, sce] = gui.gui_getfigsce(src);
    gsorted = natsort(sce.g);

    rng("shuffle");
    n = length(gsorted);
    if isempty(predefinedlist)
        if isempty(gsorted)
            ingenelist = gui.i_inputgenelist(sprintf("A1CF\nA4GNT" + ...
                "\nAADAC\nABCC8\nACADL\nACKR1\nACMSD\nACSL6\nACSM2B" + ...
                "\nACSM3\nACSM6\nACTL6B\nADAMTS8\nADCYAP1\nADGRV1\n" + ...
                "ADH1A\nADH1B\nADIPOQ\nAGR3\nAGTR2\nAJAP1\nAKAIN1\n" + ...
                "AKR7A3\nAKR7L\nALB\nALDOB\nALKAL2\nAMBP\nAMHR2\nAMPD1" + ...
                "\nAMY1B\nAMY2A\nAMY2B\nANKRD20A1\nANKRD34C\nANKRD40CL" + ...
                "\nANKRD62\nANKS4B\nANO5\nANPEP\nAPOBEC2\nAPOH\nAQP12A\nAQP12B\nAQP8\nARHGDIG\nARSL\nART3\nART4\nARX\nASB16\nASCL5\nASTN1\nATP1A2\nATP2A3\nATRNL1\nAVPR1B\nB3GAT1\nB3GNT6\nBANF2\nBEGAIN\nBEX1\nBHLHA15\nBLK\nBMP5\nBNIP5\nBRSK2\nBTBD17\nBTNL3\nBTNL8\nBTNL9\nC10orf62\nC10orf71\nC12orf42\nC14orf180\nC16orf89\nC1orf127\nC1QL1\nC22orf42\nC2CD4B\nC2orf72\nC4BPA\nC6\nC7\nCA4\nCABP7\nCACNA1B\nCACNA2D2\nCALCB\nCALN1\nCAPN6\nCAPN9\nCARTPT\nCASR\nCBFA2T3\nCBLIF\nCCDC141\nCCDC196\nCCDC92B\nCCKBR\nCCL14\nCCL15\nCCL19\nCCL21\nCCR7\nCD160\nCD27\nCD300LG\nCD8A\nCDH17\nCDH18\nCDH22\nCDHR5\nCDX2\nCEL\nCELA2A\nCELA2B\nCELA3A\nCELA3B\nCELF3\nCFAP97D2\nCFC1\nCFTR\nCHAD\nCHGA\nCHGB\nCHRDL1\nCHRNB2\nCHRNB3\nCHST5\nCHST8\nCHST9\nCIDEC\nCLCNKA\nCLDN18\nCLDN3\nCLEC17A\nCLEC9A\nCLPS\nCLRN3\nCLU\nCMA1\nCNMD\nCNR1\nCNTFR\nCNTN2\nCOLEC11\nCPA1\nCPA2\nCPB1\nCPLX2\nCR2\nCRH\nCRP\nCRYBA2\nCTNND2\nCTRB1\nCTRB2\nCTRC\nCTRL\nCTXND1\nCUX2\nCUZD1\nCXCL13\nDACH2\nDCDC2\nDCX\nDDC\nDLK1\nDMRTC1B\nDNASE1L3\nDNASE2B\nDPEP1\nDPP10\nDPYS\nDPYSL5\nDUSP26\nECE2\nEDN3\nEGF\nELAPOR1\nENAM\nENHO\nENTPD8\nEPB41L4B\nEPCIP\nEPHA8\nEPO\nERBB4\nERICH5\nERP27\nF11\nF5\nFABP4\nFAM107A\nFAM177B\nFAM180B\nFAM240C\nFAM3B\nFBP2\nFBXW12\nFCRL3\nFFAR1\nFGFR4\nFGL1\nFLT3\nFMO5\nFOXA2\nFOXI1\nFRMPD1\nFSTL5\nFUT9\nFXYD2\nFYB2\nG6PC1\nG6PC2\nGABRA4\nGABRB1\nGABRG1\nGAD2\nGALNT16\nGATM\nGBA3\nGC\nGCG\nGCGR\nGCNT3\nGFAP\nGHRL\nGHSR\nGJC3\nGJD2\nGLRA1\nGLRA3\nGLS2\nGMNC\nGNAT3\nGNMT\nGP2\nGPA33\nGPBAR1\nGPD1\nGPHA2\nGPM6A\nGPR119\nGPR142\nGPR148\nGPT2\nGRB14\nGRIA2\nGRIA4\nGRIK3\nGRPR\nGSTA2\nGSTM1\nGUCA1C\nGUCA2A\nGUCY2C\nGZMK\nHABP2\nHEPACAM2\nHHATL\nHLF\nHMGCLL1\nHMX3\nHNF1B\nHOGA1\nHOXB13\nHPD\nHTR1A\nIAPP\nIGFALS\nIGFN1\nIGLL5\nIGSF1\nIHH\nIL12B\nIL22RA1\nINS\nINS-IGF2\nINSC\nINSM1\nIRX2\nITGAD\nITIH4\nITLN1\nITPRID1\nJCHAIN\nJSRP1\nKASH5\nKCNA2\nKCNA5\nKCNC1\nKCNG3\nKCNH6\nKCNH7\nKCNIP1\nKCNJ11\nKCNJ16\nKCNJ3\nKCNK10\nKCNK16\nKCNK3\nKCNMB2\nKCNQ2\nKCNT1\nKCTD8\nKIF12\nKIF19\nKIF1A\nKIRREL2\nKLB\nKLF15\nKLHL1\nKLHL32\nKLK1\nKLKB1\nKLRB1\nKNG1\nKRT20\nKRT222\nLCN10\nLCN6\nLEFTY1\nLGALS2\nLHFPL4\nLILRA4\nLMO3\nLRRC19\nLRRC66\nLRRC7\nLRRTM3\nLRRTM4\nLTF\nLTK\nLY9\nMAFA\nMAP3K15\nMBOAT4\nMEP1A\nMFSD6L\nMGAM2\nMLXIPL\nMOGAT3\nMRLN\nMS4A1\nMT1G\nMT1H\nMTMR7\nMTUS2\nMUC13\nMUC17\nMUC5AC\nMUC6\nMUSK\nMYH7\nMYLK2\nMYMX\nMYO1A\nMYO7B\nMYOT\nMYPN\nMYRFL\nMYRIP\nMYT1\nMYT1L\nMZB1\nNAT16\nNAT8L\nNEURL1\nNEURL3\nNEUROD1\nNHERF4\nNKX2-2\nNKX6-1\nNKX6-2\nNKX6-3\nNLRP2\nNLRP9\nNOL4\nNOS1\nNPAP1\nNPHS1\nNPY\nNR0B1\nNR0B2\nNR1H4\nNR5A2\nNRAP\nNRG4\nNRSN1\nNRXN1\nNT5C1A\nNTRK2\nNXPE1\nOGDHL\nOLFM1\nOMG\nONECUT1\nOPRM1\nOR10G4\nOR10G8\nOR10G9\nOR2T11\nOR4D5\nOR52N1\nOR8D4\nORM2\nOTC\nP2RX1\nP2RX2\nPAH\nPAK3\nPAK5\nPAX5\nPCDHA1\nPCP4\nPCSK1\nPCSK1N\nPCSK2\nPDIA2\nPDK4\nPDX1\nPDZK1\nPEX5L\nPGA4\nPGC\nPHGR1\nPHYHIPL\nPIGR\nPIP5K1B\nPIRT\nPLA2G10\nPLA2G1B\nPLA2G2A\nPLA2G2D\nPLIN5\nPLPPR1\nPM20D1\nPNLIP\nPNLIPRP1\nPNLIPRP2\nPOTEC\nPOU6F2\nPPARGC1A\nPPP1R1A\nPPP1R1B\nPPY\nPRLHR\nPRODH2\nPROX1\nPRPH\nPRRT1B\nPRSS1\nPRSS2\nPRSS3\nPTF1A\nPTPN5\nPTPRN2\nPVRIG\nRAB3C\nRALYL\nRBP2\nRBP4\nRBPJL\nREG1A\nREG1B\nREG3A\nREG3G\nREG4\nRELN\nRERGL\nRFX6\nRIC3\nRIMBP2\nRIMS1\nRNF186\nRPL3L\nRXRG\nSBK2\nSCARA5\nSCG3\nSCGN\nSCML4\nSCRT1\nSERPINA10\nSERPINA6\nSERPINI2\nSEZ6\nSEZ6L\nSFRP5\nSGCZ\nSGSM1\nSH2D6\nSH2D7\nSH3GL2\nSHISAL2B\nSI\nSIM1\nSLC12A1\nSLC16A12\nSLC17A4\nSLC17A6\nSLC1A2\nSLC22A31\nSLC25A18\nSLC26A5\nSLC28A2\nSLC2A2\nSLC30A2\nSLC30A8\nSLC38A11\nSLC38A3\nSLC39A5\nSLC3A1\nSLC4A10\nSLC4A4\nSLC5A8\nSLC5A9\nSLC6A17\nSLC6A19\nSLC6A4\nSLC7A14\nSLC7A2\nSLC8A2\nSLITRK1\nSLITRK5\nSMIM31\nSMIM32\nSMLR1\nSNTG2\nSORCS1\nSOWAHA\nSPAG6\nSPATA21\nSPIB\nSPINK1\nSPINK4\nSPTB\nSPX\nSST\nSSTR1\nSSTR3\nSSTR5\nST8SIA3\nSTAB2\nSTMND1\nSULT1C2\nSV2B\nSVOP\nSYBU\nSYCN\nSYNPR\nSYT4\nSYT6\nTAGLN3\nTCEAL2\nTCERG1L\nTCL1A\nTDRD1\nTENT5C\nTEX11\nTHBS4\nTHRSP\nTM4SF20\nTM4SF5\nTMC2\nTMED6\nTMEFF2\nTMEM132C\nTMEM132D\nTMEM179\nTMEM196\nTMEM238L\nTMEM52\nTMEM63C\nTMEM72\nTNFRSF17\nTNR\nTNXB\nTPO\nTPSD1\nTRARG1\nTRIM50\nTRIM63\nTRPM3\nTRPM5\nTRPV6\nTSPAN8\nTSPEAR\nTTLL6\nTTN\nTTR\nTUBB4A\nTUNAR\nUBD\nUCN3\nUMOD\nUNC5CL\nUNC5D\nUNC79\nUNC80\nVIL1\nVIPR2\nVSTM2A\nVTN\nVWA5B2\nWNK2\nXCR1\nXKR4\nZBTB16\nZDHHC22\nZG16\nZNF831"));
        else
            ingenelist = gui.i_inputgenelist(gsorted(randperm(n, ...
                min([200, length(gsorted)]))));
        end
    else
        % ingenelist = gui.i_inputgenelist(predefinedlist);
        ingenelist = predefinedlist;
    end
    if isempty(ingenelist) || all(strlength(ingenelist) < 1)
        return;
    end


    if askbackground && n > 100
        answer = gui.myQuestdlg(FigureHandle, 'Add background list?','');
        switch answer
            case 'Yes'
                [idx] = gui.i_selmultidlg(sce.g, sce.g, FigureHandle);
                if isempty(idx), return; end
                if idx == 0, return; end
                backgroundlist = sce.g(idx);
            case 'No'
                backgroundlist = [];
            case 'Cancel'
                return;
        end
    else
        backgroundlist = [];
    end

    if isempty(enrichrtype)
        answer1 = gui.myQuestdlg(FigureHandle, ...
              "Select the type of Enrichr application.","", ...
             {'Web-based', 'API-based'}, 'API-based');
        enrichrtype = answer1;
    end

switch enrichrtype 
    case 'API-based'
        % do nothing here


    case 'Web-based'
        fw = gui.myWaitbar(FigureHandle, [], false, ...
            'Sending genes to web browser...');
        % gui.i_enrichtest(genelist, backgroundlist, numel(genelist));
            if ~isempty(backgroundlist)
                run.web_Enrichr_bkg(ingenelist, backgroundlist, numel(ingenelist), wkdir);
            else
                run.web_Enrichr(ingenelist, numel(ingenelist), '', wkdir);
            end
        gui.myWaitbar(FigureHandle, fw, false, ...
            'Check web browser & submit genes to Enrichr.');
        return;   
    otherwise
        return;
end



    [genesets] = in_selDataSources;
    if isempty(genesets), return; end

    definput = {'5', '0.1'};
    prompt = {'Min # of overlapping genes:', ...
              'P-value cutoff:'};
    dlgtitle = 'Enrichr Result Filter';
    dims = [1, 80];
    if gui.i_isuifig(FigureHandle)
        answer = gui.myInputdlg(prompt, dlgtitle, definput, FigureHandle);
    else
        answer = inputdlg(prompt, dlgtitle, dims, definput);
    end

    if isempty(answer), return; end
    try
        minugenes = str2double(answer{1});
        pvaluecut = str2double(answer{2});
        assert((minugenes > 0) && (minugenes < 100));
        assert((pvaluecut >= 0.0) && (pvaluecut <= 1.0));
    catch
        gui.myErrordlg(FigureHandle, 'Invalid input.');
        return;
    end
    

    fw = gui.myWaitbar(FigureHandle);
    Tlist = run.ml_Enrichr(ingenelist, backgroundlist, genesets,...
                           minugenes, pvaluecut);

    T=table;
    for k = 1:height(Tlist)
        if ~isempty(Tlist{k}) && istable(Tlist{k})
            T = [T; Tlist{k}];
        end
    end

    gui.myWaitbar(FigureHandle, fw);
    
    %[~, ~] = gui.i_exporttable(T, true, 'Tenrichrres', ...
    %    sprintf('Enrichr_Results_%s', outfiletag),[],[],FigureHandle);

    options = {'View Table', 'Circos Plot'};
    %answer = gui.myQuestdlg(FigureHandle, 'View Enrichr Result Table or Show the Table as a Circos Plot?', ...
    %    '', {options{1}, options{2}}, options{1});

    answer = options{1};
    switch answer
        case options{1}
            % gui.i_viewtable(T, FigureHandle);
            gui.TableViewerApp(T, FigureHandle);

        case options{2}
            if height(T) > 1
                gui.callback_EnrichrTab2Circos(src, [], T);
            else
                gui.myHelpdlg(FigureHandle, 'Too few enriched terms to generate circos plot.','');                
            end
        otherwise

    end


    function [genesets] = in_selDataSources
        genesets = [];
        dsv = pkg.i_get_enrichr_libraries;

        enrichrlibraries = getpref('scgeatoolbox', 'enrichrlibraries', ...
                                   ["GO_Biological_Process_2025", ...
                                    "GO_Molecular_Function_2025", ...
                                     "KEGG_2021_Human",...
                                     "Reactome_Pathways_2024"]);

        [idx1] = gui.i_selmultidlg(dsv, enrichrlibraries, FigureHandle);
        if isempty(idx1), return; end
        if idx1 == 0, return; end
        genesets = dsv(idx1);
        setpref('scgeatoolbox', 'enrichrlibraries', genesets);    
    end

end
