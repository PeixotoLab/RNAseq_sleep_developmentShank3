function EEG_struct = SetUpEEGstruct(EEGWTBLWake,   EEGWTBLNREM,   EEGWTBLREM,   EEGWTSDWake,   EEGWTSDNREM,   EEGWTSDREM,...
                            		 EEGMutBLWake,  EEGMutBLNREM,  EEGMutBLREM,  EEGMutSDWake,  EEGMutSDNREM,  EEGMutSDREM, ...
                            		 EEGWTBLWake_F, EEGWTBLNREM_F, EEGWTBLREM_F, EEGWTSDWake_F, EEGWTSDNREM_F, EEGWTSDREM_F,...
                            		 EEGMutBLWake_F,EEGMutBLNREM_F,EEGMutBLREM_F,EEGMutSDWake_F,EEGMutSDNREM_F,EEGMutSDREM_F,Sexes)

% This helper function sets up the struct EEG_struct. It cleans up the main script just a bit.  



% - Spectral Data -
% Males  (Spectral Data)  
EEG_struct.WT.(Sexes{1}).BL.Wake = EEGWTBLWake;  % WT
EEG_struct.WT.(Sexes{1}).BL.NREM = EEGWTBLNREM;
EEG_struct.WT.(Sexes{1}).BL.REM  = EEGWTBLREM;
EEG_struct.WT.(Sexes{1}).SD.Wake = EEGWTSDWake;
EEG_struct.WT.(Sexes{1}).SD.NREM = EEGWTSDNREM;
EEG_struct.WT.(Sexes{1}).SD.REM  = EEGWTSDREM;

EEG_struct.Mut.(Sexes{1}).BL.Wake = EEGMutBLWake;  % Mutant
EEG_struct.Mut.(Sexes{1}).BL.NREM = EEGMutBLNREM;
EEG_struct.Mut.(Sexes{1}).BL.REM  = EEGMutBLREM;
EEG_struct.Mut.(Sexes{1}).SD.Wake = EEGMutSDWake;
EEG_struct.Mut.(Sexes{1}).SD.NREM = EEGMutSDNREM;
EEG_struct.Mut.(Sexes{1}).SD.REM  = EEGMutSDREM;

% Females  (Spectral Data)  
EEG_struct.WT.(Sexes{2}).BL.Wake = EEGWTBLWake_F;  % WT
EEG_struct.WT.(Sexes{2}).BL.NREM = EEGWTBLNREM_F;
EEG_struct.WT.(Sexes{2}).BL.REM  = EEGWTBLREM_F;
EEG_struct.WT.(Sexes{2}).SD.Wake = EEGWTSDWake_F;
EEG_struct.WT.(Sexes{2}).SD.NREM = EEGWTSDNREM_F;
EEG_struct.WT.(Sexes{2}).SD.REM  = EEGWTSDREM_F;

EEG_struct.Mut.(Sexes{2}).BL.Wake = EEGMutBLWake_F;  % Mutant
EEG_struct.Mut.(Sexes{2}).BL.NREM = EEGMutBLNREM_F;
EEG_struct.Mut.(Sexes{2}).BL.REM  = EEGMutBLREM_F;
EEG_struct.Mut.(Sexes{2}).SD.Wake = EEGMutSDWake_F;
EEG_struct.Mut.(Sexes{2}).SD.NREM = EEGMutSDNREM_F;
EEG_struct.Mut.(Sexes{2}).SD.REM  = EEGMutSDREM_F;
