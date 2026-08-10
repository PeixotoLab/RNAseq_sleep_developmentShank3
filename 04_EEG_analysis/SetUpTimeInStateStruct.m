function TimeInState = SetUpTimeInStateStruct(Wake24hWTBL,   NREM24hWTBL,   REM24hWTBL,   Wake24hWTSD,   NREM24hWTSD,   REM24hWTSD, ...
                                     		  Wake24hMutBL,  NREM24hMutBL,  REM24hMutBL,  Wake24hMutSD,  NREM24hMutSD,  REM24hMutSD, ...
                                     		  Wake24hWTBL_F, NREM24hWTBL_F, REM24hWTBL_F, Wake24hWTSD_F, NREM24hWTSD_F, REM24hWTSD_F, ...
                                     		  Wake24hMutBL_F,NREM24hMutBL_F,REM24hMutBL_F,Wake24hMutSD_F,NREM24hMutSD_F,REM24hMutSD_F,SD_length_hrs,Sexes)

% This is just a little helper function to set up the TimeInState struct so the main function isn't too cluttered. 





% WT Male 
TimeInState.WT.(Sexes{1}).BL.Wake = Wake24hWTBL;
TimeInState.WT.(Sexes{1}).BL.NREM = NREM24hWTBL;
TimeInState.WT.(Sexes{1}).BL.REM  = REM24hWTBL;
TimeInState.WT.(Sexes{1}).SD.Wake = Wake24hWTSD;                                % SD means the day including the sleep dep (the entire day)
TimeInState.WT.(Sexes{1}).SD.NREM = NREM24hWTSD;      
TimeInState.WT.(Sexes{1}).SD.REM  = REM24hWTSD; 
TimeInState.WT.(Sexes{1}).SDexcSD.Wake = Wake24hWTSD(SD_length_hrs+1:end,:);    % this is the sleep dep day with the actual sleep dep excluded
TimeInState.WT.(Sexes{1}).SDexcSD.NREM = NREM24hWTSD(SD_length_hrs+1:end,:);
TimeInState.WT.(Sexes{1}).SDexcSD.REM  = REM24hWTSD(SD_length_hrs+1:end,:);

% Mutant Male 
TimeInState.Mut.(Sexes{1}).BL.Wake = Wake24hMutBL;
TimeInState.Mut.(Sexes{1}).BL.NREM = NREM24hMutBL;
TimeInState.Mut.(Sexes{1}).BL.REM  = REM24hMutBL;
TimeInState.Mut.(Sexes{1}).SD.Wake = Wake24hMutSD;                                % SD means the day including the sleep dep (the entire day)
TimeInState.Mut.(Sexes{1}).SD.NREM = NREM24hMutSD;      
TimeInState.Mut.(Sexes{1}).SD.REM  = REM24hMutSD; 
TimeInState.Mut.(Sexes{1}).SDexcSD.Wake = Wake24hMutSD(SD_length_hrs+1:end,:);    % this is the sleep dep day with the actual sleep dep excluded
TimeInState.Mut.(Sexes{1}).SDexcSD.NREM = NREM24hMutSD(SD_length_hrs+1:end,:);
TimeInState.Mut.(Sexes{1}).SDexcSD.REM  = REM24hMutSD(SD_length_hrs+1:end,:);

% WT Female 
TimeInState.WT.(Sexes{2}).BL.Wake = Wake24hWTBL_F;
TimeInState.WT.(Sexes{2}).BL.NREM = NREM24hWTBL_F;
TimeInState.WT.(Sexes{2}).BL.REM  = REM24hWTBL_F;
TimeInState.WT.(Sexes{2}).SD.Wake = Wake24hWTSD_F;                                % SD means the day including the sleep dep (the entire day)
TimeInState.WT.(Sexes{2}).SD.NREM = NREM24hWTSD_F;      
TimeInState.WT.(Sexes{2}).SD.REM  = REM24hWTSD_F; 
TimeInState.WT.(Sexes{2}).SDexcSD.Wake = Wake24hWTSD_F(SD_length_hrs+1:end,:);    % this is the sleep dep day with the actual sleep dep excluded
TimeInState.WT.(Sexes{2}).SDexcSD.NREM = NREM24hWTSD_F(SD_length_hrs+1:end,:);
TimeInState.WT.(Sexes{2}).SDexcSD.REM  = REM24hWTSD_F(SD_length_hrs+1:end,:);

% Mutant Female 
TimeInState.Mut.(Sexes{2}).BL.Wake = Wake24hMutBL_F;
TimeInState.Mut.(Sexes{2}).BL.NREM = NREM24hMutBL_F;
TimeInState.Mut.(Sexes{2}).BL.REM  = REM24hMutBL_F;
TimeInState.Mut.(Sexes{2}).SD.Wake = Wake24hMutSD_F;                                % SD means the day including the sleep dep (the entire day)
TimeInState.Mut.(Sexes{2}).SD.NREM = NREM24hMutSD_F;      
TimeInState.Mut.(Sexes{2}).SD.REM  = REM24hMutSD_F; 
TimeInState.Mut.(Sexes{2}).SDexcSD.Wake = Wake24hMutSD_F(SD_length_hrs+1:end,:);    % this is the sleep dep day with the actual sleep dep excluded
TimeInState.Mut.(Sexes{2}).SDexcSD.NREM = NREM24hMutSD_F(SD_length_hrs+1:end,:);
TimeInState.Mut.(Sexes{2}).SDexcSD.REM  = REM24hMutSD_F(SD_length_hrs+1:end,:);

