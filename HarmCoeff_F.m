clear all; close all; clc
mu0=4*pi*1e-7;
nSP=12;
p=100;
q=40;
mesh=load('J_100by40 per SP.txt');
mesh_sorted=sortrows(round(mesh.*1e8)./1e8,[2,1]);
mesh_R=reshape(mesh_sorted(:,1).*1e-3,[(p+1),(q+1),nSP]);
mesh_Z=reshape(mesh_sorted(:,2).*1e-3,[(p+1),(q+1),nSP]);
% mesh_J=reshape(mesh_sorted(:,3),[(p+1),(q+1),nSP]);
mesh_J=reshape(mesh_sorted(:,3),[(p+1),(q+1)*nSP]);
mesh_J=(mesh_J+fliplr(mesh_J))./2;
mesh_J=reshape(mesh_J,[(p+1),(q+1),nSP]);
%%
Area1=zeros(p,q,nSP);
cenR1=zeros(p,q,nSP);
cenZ1=zeros(p,q,nSP);
cenJ1=zeros(p,q,nSP);
Area2=zeros(p,q,nSP);
cenR2=zeros(p,q,nSP);
cenZ2=zeros(p,q,nSP);
cenJ2=zeros(p,q,nSP);
for j=1:nSP
    for k=1:q
        for l=1:p
            %         |-/
            x1=mesh_R(l,k,j);        y1=mesh_Z(l,k,j);        w1=mesh_J(l,k,j);
            x2=mesh_R(l+1,k,j);      y2=mesh_Z(l+1,k,j);      w2=mesh_J(l+1,k,j);
            x3=mesh_R(l+1,k+1,j);    y3=mesh_Z(l+1,k+1,j);    w3=mesh_J(l+1,k+1,j);
            Area1(l,k,j)=ShoeLaceFormula(x1,y1,x2,y2,x3,y3);
            [xg,yg,wg]=CenterMassTriangle(x1,y1,w1,x2,y2,w2,x3,y3,w3);
            cenR1(l,k,j)=xg;        cenZ1(l,k,j)=yg;        cenJ1(l,k,j)=wg;
            %         /_|
            x1=mesh_R(l,k,j);        y1=mesh_Z(l,k,j);        w1=mesh_J(l,k,j);
            x2=mesh_R(l,k+1,j);      y2=mesh_Z(l,k+1,j);      w2=mesh_J(l,k+1,j);
            x3=mesh_R(l+1,k+1,j);    y3=mesh_Z(l+1,k+1,j);    w3=mesh_J(l+1,k+1,j);
            Area1(l,k,j)=Area1(l,k,j)+1i*ShoeLaceFormula(x1,y1,x2,y2,x3,y3);
            [xg,yg,wg]=CenterMassTriangle(x1,y1,w1,x2,y2,w2,x3,y3,w3);
            cenR1(l,k,j)=cenR1(l,k,j)+1i*xg;        cenZ1(l,k,j)=cenZ1(l,k,j)+1i*yg;        cenJ1(l,k,j)=cenJ1(l,k,j)+1i*wg;
            %         \-|
            x1=mesh_R(l,k,j);        y1=mesh_Z(l,k,j);        w1=mesh_J(l,k,j);
            x2=mesh_R(l,k+1,j);      y2=mesh_Z(l,k+1,j);      w2=mesh_J(l,k+1,j);
            x3=mesh_R(l+1,k,j);      y3=mesh_Z(l+1,k,j);      w3=mesh_J(l+1,k,j);
            Area2(l,k,j)=ShoeLaceFormula(x1,y1,x2,y2,x3,y3);
            [xg,yg,wg]=CenterMassTriangle(x1,y1,w1,x2,y2,w2,x3,y3,w3);
            cenR2(l,k,j)=xg;        cenZ2(l,k,j)=yg;        cenJ2(l,k,j)=wg;
            %         |_\
            x1=mesh_R(l+1,k,j);        y1=mesh_Z(l+1,k,j);        w1=mesh_J(l+1,k,j);
            x2=mesh_R(l,k+1,j);        y2=mesh_Z(l,k+1,j);        w2=mesh_J(l,k+1,j);
            x3=mesh_R(l+1,k+1,j);      y3=mesh_Z(l+1,k+1,j);      w3=mesh_J(l+1,k+1,j);
            Area2(l,k,j)=Area2(l,k,j)+1i*ShoeLaceFormula(x1,y1,x2,y2,x3,y3);
            [xg,yg,wg]=CenterMassTriangle(x1,y1,w1,x2,y2,w2,x3,y3,w3);
            cenR2(l,k,j)=cenR2(l,k,j)+1i*xg;        cenZ2(l,k,j)=cenZ2(l,k,j)+1i*yg;        cenJ2(l,k,j)=cenJ2(l,k,j)+1i*wg;
        end
    end
end
Area1=reshape(Area1,[p,q*nSP]);
cenR1=reshape(cenR1,[p,q*nSP]);
cenZ1=reshape(cenZ1,[p,q*nSP]);
cenJ1=reshape(cenJ1,[p,q*nSP]);
Area2=reshape(Area2,[p,q*nSP]);
cenR2=reshape(cenR2,[p,q*nSP]);
cenZ2=reshape(cenZ2,[p,q*nSP]);
cenJ2=reshape(cenJ2,[p,q*nSP]);
cenI1=(real(Area1).*real(cenJ1))+1i*(imag(Area1).*imag(cenJ1));
cenI2=(real(Area2).*real(cenJ2))+1i*(imag(Area2).*imag(cenJ2));

% figure;
% subplot(1,2,1);
% surf(real(cenR2),real(cenZ2),real(Area2));
% axis tight; axis equal; view(0,90);
% colormap jet; colorbar;
% shading interp;
% subplot(1,2,2);
% surf(imag(cenR2),imag(cenZ2),imag(Area2)); 
% axis tight; axis equal; view(0,90);
% colormap jet; colorbar;
% shading interp;

% Harmonic Coefficients w. Associate Legendre
nOrder=10;
ALegen1=zeros(p, q*nSP ,nOrder+1);
for j=1:p
    for k=1:q*nSP
        Isrc1=real(cenI1(j,k));
        Rsrc1=real(cenR1(j,k));
        Zsrc1=real(cenZ1(j,k));
        c_val1=(Zsrc1)/sqrt(Rsrc1^2+Zsrc1^2);
        s_val1=(Rsrc1)/sqrt(Rsrc1^2+Zsrc1^2);
        Isrc2=imag(cenI1(j,k));
        Rsrc2=imag(cenR1(j,k));
        Zsrc2=imag(cenZ1(j,k));
        c_val2=(Zsrc2)/sqrt(Rsrc2^2+Zsrc2^2);
        s_val2=(Rsrc2)/sqrt(Rsrc2^2+Zsrc2^2);
        for l=0:nOrder
            ALegen1(j,k,l+1)=-1*AssoLegen_Coeff(1,0,l)*(mu0/2)*Isrc1*(s_val1/sqrt(Rsrc1^2+Zsrc1^2)^(l+1))*AssoLegen_Coeff(c_val1,1,l+1)+...
                -1*AssoLegen_Coeff(1,0,l)*(mu0/2)*Isrc2*(s_val2/sqrt(Rsrc2^2+Zsrc2^2)^(l+1))*AssoLegen_Coeff(c_val2,1,l+1);
        end
    end
end

ALegen2=zeros(p, q*nSP,nOrder+1);
for j=1:p
    for k=1:q*nSP
        Isrc1=real(cenI2(j,k));
        Rsrc1=real(cenR2(j,k));
        Zsrc1=real(cenZ2(j,k));
        c_val1=(Zsrc1)/sqrt(Rsrc1^2+Zsrc1^2);
        s_val1=(Rsrc1)/sqrt(Rsrc1^2+Zsrc1^2);
        Isrc2=imag(cenI2(j,k));
        Rsrc2=imag(cenR2(j,k));
        Zsrc2=imag(cenZ2(j,k));
        c_val2=(Zsrc2)/sqrt(Rsrc2^2+Zsrc2^2);
        s_val2=(Rsrc2)/sqrt(Rsrc2^2+Zsrc2^2);
        for l=0:nOrder
            ALegen2(j,k,l+1)=-1*AssoLegen_Coeff(1,0,l)*(mu0/2)*Isrc1*(s_val1/sqrt(Rsrc1^2+Zsrc1^2)^(l+1))*AssoLegen_Coeff(c_val1,1,l+1)+...
                -1*AssoLegen_Coeff(1,0,l)*(mu0/2)*Isrc2*(s_val2/sqrt(Rsrc2^2+Zsrc2^2)^(l+1))*AssoLegen_Coeff(c_val2,1,l+1);
        end
    end
end

%%
HarmCoeff=zeros(nOrder+1,1);
HarmCoeff_cm=zeros(nOrder+1,1);
HarmCoeff_mm=zeros(nOrder+1,1);
Zn=zeros(nOrder+1,1);
Zn_cm=zeros(nOrder+1,1);
Zn_mm=zeros(nOrder+1,1);
for j=1:nOrder+1
    HarmCoeff(j)=sum(sum(squeeze(ALegen1(1:end,1:end,j)+ALegen2(1:end,1:end,j))))./2;
    HarmCoeff_cm(j)=HarmCoeff(j)/(100^(j-1));
    HarmCoeff_mm(j)=HarmCoeff(j)/(1000^(j-1));
    Zn(j)=HarmCoeff(j)*factorial(j-1);
    Zn_cm(j)=Zn(j)/100^(j-1);
    Zn_mm(j)=Zn(j)/1000^(j-1);
end

data=load('Bz_100by40 per SP_SCF.txt');
z1=data(:,2); z2=-data(:,2);
Bz=data(:,3);
zz=linspace(-10,10,101);
Bzz=(interp1(z1,Bz,zz)+interp1(z2,Bz,zz))./2;
a_mm=[0.95;0;-0.000286;0;4.213e-8;0;-1.463e-11;0;5.193e-14;0;-5.511e-17];
a_cm=Zn_cm./HarmCoeff_mm.*a_mm;