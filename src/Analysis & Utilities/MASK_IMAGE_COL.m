function [I_MASKED] = MASK_IMAGE_COL(I,M,COL)

[h,w] = size(I);

TEMP1 = I; 
TEMP2 = I; 
TEMP3 = I;
TEMP1(M==0) = COL(1); 
TEMP2(M==0) = COL(2); 
TEMP3(M==0) = COL(3);
I_MASKED        = zeros(h,w,3); 
I_MASKED(:,:,1) = TEMP1; 
I_MASKED(:,:,2) = TEMP2; 
I_MASKED(:,:,3) = TEMP3; 


end

