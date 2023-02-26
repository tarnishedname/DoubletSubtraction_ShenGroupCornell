#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.



//Hello! Here's an example of how to use this code: 
//"root:Example_NeI_Data_Bi2Se3:NeI_presubtracted" is an example ARPES images of single crystal Bi2Se3
//Note that is has been converted from angle to momentum on the x-axis.
//Set "Example_NeI_Data_Bi2Se3" as your current data folder
//Run "BetaSubtraction_(NeI_presubtracted, 0.176, 0.25)"


function BetaSubtraction_(wave image, variable doublet_dE, variable doublet_ratio)
//Input ARPES image with x-axis as k, y-axis as E-EF. Input the accurately-determined dE. Input a reasonable guess to R

variable i = 0
Make/O/N=(dimsize(image,1)) IntegratedEDC
Setscale/P x,dimoffset(image,1), dimdelta(image,1), IntegratedEDC

for(i=0;i<dimsize(image,1);i++)
  Duplicate/O/R=[][i]/FREE image image_temp_
  Wavestats/Q image_temp_
  IntegratedEDC[i] = V_avg 
endfor   

Duplicate/O/R=(0.05,inf) IntegratedEDC autoLinBG  //This part subtracts a linear background starting slightly above the Fermi level.
CurveFit line autoLinBG /D 
wave fit_autoLinBG
wave W_coef
image-=(W_coef[0]+W_coef[1]*y)


variable optimal_ratio = Optimize_R_no_disc(IntegratedEDC, doublet_dE, doublet_ratio) //Calls function to determine the optimal R value.
variable full_k_min = dimOffset(image,0)
variable full_k_max = dimOffset(image,0)+dimDelta(image,0)*(dimSize(image,0)-1)
Make/O/N=(200,1) FullImageFitted //Creates space for output image. Size doesn't matter, it's redimensioned later anyway.
variable EDCindex = 0
i=0		   
for(i=full_k_min;i<=(full_k_max);i+=((full_k_max-full_k_min)/dimsize(image,0))) //Sampling to same number as original.
	EDCindex+=1
	wave full_extract = auto_ExtractEDCs(image,i)
	Main_Doublet_Algorithm(full_extract, 0, doublet_dE, optimal_ratio) //Calls main algorithm, feeding it the optimal R determined above.
	wave returned_noff_fit
	Redimension/N=(dimsize(image,1),EDCindex) FullImageFitted
	Setscale/P x,dimoffset(image,1), dimdelta(image,1), FullImageFitted
	SetScale/I y, dimOffset(image,0),dimOffset(image,0)+dimDelta(image,0)*(dimSize(image,0)-1), FullImageFitted
	Duplicate/O returned_noff_fit fullfitsappend
	FullImageFitted[][EDCindex-1] = fullfitsappend[p] //Appends all the post-subtracted EDCs together to create output image		       
endfor		      		   
Duplicate/O image image_transpose_local
MatrixTranspose image_transpose_local
auto_compare_images(image_transpose_local,FullImageFitted) //For display and diagnostic purposes
print(optimal_ratio)

KillWaves/Z fit_sub_returned fullfitsappend image_transpose_local Res_sub_returned returned_noff_fit sub_returned test_local test_local_sub test_local_subshift  

end



function Optimize_R_no_disc(wave AverageEDC, variable dE, variable ratio)
//This function takes the integrated EDC of the spectra, the accurately-determined DeltaE, and a fair guess of R, and outputs the optimal R
variable i
variable ratio_0
variable foundBool=0
Make/O/D/N=(411,2) chi_sq_vs_ratio
for(i=0;i<411;i+=1)
  ratio_0 = (ratio-0.1)+i*(0.0005)
  Main_Doublet_Algorithm(AverageEDC, 0, dE, ratio_0)
  wave returned_noff_fit
  Duplicate/O/R=(-1.5*dE,-0.5*dE) returned_noff_fit sub_returned  //Range of either +-25% OR +-50% seems appropriate.
  CurveFit/Q poly 3, sub_returned /D /R /A=1  //Smooth polynomial fit to straddle the discontinuity.
  variable chi_sq = V_chisq
  chi_sq_vs_ratio[i][0] = ratio_0
  chi_sq_vs_ratio[i][1] = chi_sq
endfor
MDsort(chi_sq_vs_ratio,1) //Sorts such that mininum chi_sq is 1st row
//display chi_sq_vs_ratio[][1] vs chi_sq_vs_ratio[][0]  //Uncomment if you want to see the chi_square vs R curve.
variable optimal_ratio = chi_sq_vs_ratio[0][0]
return optimal_ratio   
End  



Function/WAVE auto_ExtractEDCs(wave2d, k_pos) //Simple function that extracts EDCs at some k. 
wave wave2d
variable k_pos
Make/O/N=(dimSize(wave2d,1))/FREE wave1d  
SetScale/P x,dimOffset(wave2d,1), dimDelta(wave2d,1), wave1d
wave1d = wave2d(k_pos)[p]

return wave1d
end


function/wave Main_Doublet_Algorithm(wave test, variable Ef, variable dE, variable ratio) //Main subtraction algorithm. See main text and Figure 2.
variable i=0
Duplicate/O test test_local

do 
   Duplicate/O/R=(Ef-dE/2-i*dE,Ef+dE/2) test_local test_local_sub 
   Duplicate/O test_local_sub test_local_subshift
   test_local_subshift*=ratio
   SetScale/P x,(dimoffset(test_local_subshift,0)-dE),dimdelta(test_local_subshift,0), test_local_subshift
   test_local = (x>dimoffset(test_local_subshift,0)) && (x<dimoffset(test_local_subshift,0)+dimdelta(test_local_subshift,0)*(dimsize(test_local_subshift,0)-1)) ? test(x)-test_local_subshift(x): test(x) 
   i+=1
  //print(i)
while (dimoffset(test_local_subshift,0)>dimoffset(test_local,0))

Duplicate/O test_local returned_noff_fit
return returned_noff_fit
End



function auto_compare_images(wave raw_minus_bg, wave fullyfitted)  //Just a function to display and compare pre and post-subtracted spectra.
Duplicate/O raw_minus_bg FullMinusBG_tau
Duplicate/O fullyfitted FullImageFitted_tau
MatrixTranspose FullMinusBG_tau
MatrixTranspose FullImageFitted_tau
NewImage/F/K=0/N=Raw_Minus_BG FullMinusBG_tau
NewImage/F/K=0/N=FullyFitted FullImageFitted_tau
ModifyImage/W=Raw_Minus_BG FullMinusBG_tau ctab={*,*,BlackBody,0} 
ModifyImage/W=FullyFitted FullImageFitted_tau ctab={*,*,BlackBody,0} //Just an arbitrary choice of color scale
ModifyGraph/W=FullyFitted nticks(bottom)=3;DelayUpdate
Label/W=FullyFitted bottom "Momentum (1/Å)"
Label/W=FullyFitted left "Energy (eV)"
ModifyGraph/W=FullyFitted margin(left)=28,margin(bottom)=28
end


Function MDsort(w,keycol, [reversed]) //Sorting function
    Wave w
    variable keycol, reversed
   
    variable type
   
    type = Wavetype(w)
 
    make/Y=(type)/free/n=(dimsize(w,0)) key
    make/free/n=(dimsize(w,0)) valindex
   
    if(type == 0)
        Wave/t indirectSource = w
        Wave/t output = key
        output[] = indirectSource[p][keycol]
    else
        Wave indirectSource2 = w
        multithread key[] = indirectSource2[p][keycol]
    endif
   
    valindex=p
    if(reversed)
        sort/a/r key,key,valindex
    else
        sort/a key,key,valindex
    endif
   
    if(type == 0)
        duplicate/free indirectSource, M_newtoInsert
        Wave/t output = M_newtoInsert
        output[][] = indirectSource[valindex[p]][q]
        indirectSource = output
    else
        duplicate/free indirectSource2, M_newtoInsert
        multithread M_newtoinsert[][] = indirectSource2[valindex[p]][q]
        multithread indirectSource2 = M_newtoinsert
    endif
End