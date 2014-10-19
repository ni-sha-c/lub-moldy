clear all
close all
%lattice constant
a0= 2.855324e0;
%Number of lattices in x and y direction
N = 10;
%System height in lattice units
H = 15;
%No of hexadecane molecules
nhd = 25;
%No of hexadecane atoms
nahd = 1250;

%Bottom Half
%Base Lattice
%Layer 1
x_latt1 = [0.e0, a0, a0, 0.e0, 0.e0, a0, a0, 0.e0, a0/2.e0];
y_latt1 = [0.e0, 0.e0, a0, a0, a0, a0, 0.e0, 0.e0, a0/2.e0];
z_latt1 = [0.e0, 0.e0, 0.e0, 0.e0, a0, a0, a0, a0, a0/2.e0];

xp = x_latt1([2,3,6,7,9]);
yp = y_latt1([2,3,6,7,9]);
zp = z_latt1([2,3,6,7,9]);

scal1 =  1:N-1;
scal  = repmat(scal1, 5, 1);
scal  = reshape(scal, 1, 5*(N-1)); 

scal = scal.*a0;

x_lattn = repmat(xp, 1, N - 1);
y_lattn = repmat(yp, 1, N - 1);
z_lattn = repmat(zp, 1, N - 1);

x_lattn = x_lattn + scal;
x1 = [x_latt1, x_lattn];
y1 = [y_latt1, y_lattn];
z1 = [z_latt1, z_lattn];

%Layer 2

x2_latt1 = x_latt1([3,4,5,6,9]);
y2_latt1 = y_latt1([3,4,5,6,9]) + a0 ;
z2_latt1 = z_latt1([3,4,5,6,9]);

scal  = repmat(scal1, 3, 1);
scal  = reshape(scal, 1, 3*(N-1));
scal  = scal.*a0;

xp2 = x2_latt1([1,4,5]);
yp2 = y2_latt1([1,4,5]);
zp2 = z2_latt1([1,4,5]);

x2_lattn = repmat(xp2, 1, N - 1);
y2_lattn = repmat(yp2, 1, N - 1);
z2_lattn = repmat(zp2, 1, N - 1);

x2_lattn = x2_lattn + scal;
x2 = [x2_latt1, x2_lattn];
y2 = [y2_latt1, y2_lattn];
z2 = [z2_latt1, z2_lattn];

%Layers 3 to N

xn = repmat(x2, 1, N-2);
yn = repmat(y2, 1, N-2);
zn = repmat(z2, 1, N-2);

scal1 = 1:N-2;
scal  = repmat(scal1, 3*(N-1) + 5, 1);
scal  = reshape(scal, 1, (N-2)*(3*(N-1)+5));
scal  = scal.*a0;

yn = yn + scal;

%Putting all the layers together

x_base = [x1, x2, xn];
y_base = [y1, y2, yn];
z_base = [z1, z2, zn];

%Number of atoms in the base
noab = length(x_base);

%raising the base to a certain height
%Total system is a cuboid of height H lattice units. 
%The leaving some units empty for lubricant.

bh = 5;
xh = repmat(x_base, 1, bh); 
yh = repmat(y_base, 1, bh);
zh = repmat(z_base, 1, bh);

scal1 = 0:bh-1;
scal  = repmat(scal1, noab, 1);
scal  = reshape(scal, 1, noab*bh);
scal  = scal.*a0;

zh = zh + scal;
%a group of lattice units is raised to the same random height 
%in the hope of avoiding needles. - version 2
rarr = rand(1,N^2);

rh1 = floor((H-2*bh)/2*rarr);

%Top layer before increasing height randomly
z_base = z_base + (bh-1)*a0;
rh = [];
for i = 1:N^2
    if i == 1
        rh = [rh, repmat(rh1(i), 1, 9)]; %9 for BCC
    elseif mod(i,N) == 1 || i <= N
        rh = [ rh, repmat(rh1(i), 1, 5)]; %starting lattice of layer
    else
        rh = [ rh, repmat(rh1(i), 1, 3)];
    end
    
end 

for i = 1:noab
    h = rh(i); 
    scal = (1:h)*a0;
    xh = [xh, repmat(x_base(i), 1, h)];
    yh = [yh, repmat(y_base(i), 1, h)];
    zh = [zh, scal + z_base(i)]; 
    
end

%Total number of atoms(including top half)
tn = 2*length(xh);

%Printing in LAMMPS format for "full" atoms
%Rigid wall is one layer. Thermostat layers 2. 
attype = repmat([ones(1,noab), 2*ones(1, 2*noab), ...
    3*ones(1, length(xh) - 3*noab)],1,2);
A = [(nahd + 1):(tn + nahd); (nhd + 1):(tn + 25) ; attype; ... 
zeros(1, tn); xh, xh; yh, yh; zh, H*a0 - zh];
fid = fopen('solid_walls.data', 'w');
fprintf(fid, '%d %d %d %f %15.8f %15.8f %15.8f\n', A);
fclose(fid);









