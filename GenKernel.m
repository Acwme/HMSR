function [res] = GenKernel(img) %img is a 4x4 patch, for calc gradient
% img = imread('zz.png');
% img = im2double(img)% *255
%1 cal gradient
tdx = [99,99,99;99,99,99;99,99,99];
tdy = [99,99,99;99,99,99;99,99,99];
for i = 1:size(img,1)-1
    for j = 1:size(img,2)-1
        tdx(i,j) = img(i,j+1)-img(i,j); 
        tdy(i,j) = img(i+1,j)-img(i,j); 
    end
end

%2 cal Ix2, IxIy, Iy2, sigma
d = 1e-3;
a = d; b = d; c = d;
for i = 1:size(tdx,1)
    for j = 1:size(tdx,2)
        a = a + tdx(i,j)*tdx(i,j);
        b = b + tdx(i,j)*tdy(i,j);
        c = c + tdy(i,j)*tdy(i,j);
    end
end
sigma = [a,b;b,c];

%3 calc lamda1, lamda2, k1, k2
lamda1 = 0.5*(a+c+sqrt((a-c)*(a-c)+4*b*b)); %????
lamda2 = 0.5*(a+c-sqrt((a-c)*(a-c)+4*b*b));
lamda2 = max(lamda2, 1e-4);
e1k = (b/(lamda1-a) + (lamda1-c)/b)/2;
e2k = (b/(lamda2-a) + (lamda2-c)/b)/2;
e1 = [e1k*sqrt(1/(e1k*e1k+1)),sqrt(1/(e1k*e1k+1))];
e2 = [e2k*sqrt(1/(e2k*e2k+1)),sqrt(1/(e2k*e2k+1))];
D_tr = 0.01; D_th = 0.005;
k_detail = 0.29; k_denoise = 4.0; K_shrink = 2.0; k_streach = 4.0;
A = 1 + sqrt((lamda1 ) / (lamda2));
D = min(max(1.0 - sqrt(lamda1) / D_tr + D_th, 0.0), 1.0);
k_1_hat = k_detail * k_streach * A;
k_2_hat = k_detail / K_shrink / A;
k1 = ((1 - D)*k_1_hat + D * k_detail*k_denoise);
k2 = ((1 - D)*k_2_hat + D * k_detail*k_denoise);
k1 = k1*k1;
k2 = k2*k2;
sigmak_ = [e1' e2'] * [k1,0;0,k2]*[e1;e2];
% sigmak_ = [e1' e2'] * [1/k1,0;0,1/k2]*[e1;e2];
sigma_ = sigmak_;
% sigma_ = sigma;%^-1;

%4 calc kernel weights
res = [99,99,99;99,99,99;99,99,99];
for i = 1:size(tdx,2)
    for j = 1:size(tdx,1)
        vec = [i-2,j-2];
        res(i,j) = exp(-0.5*(vec*sigma_*vec'));
    end
end
end


