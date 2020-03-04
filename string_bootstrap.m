n = 100; %Order of taylor expansion
D = 4; %Dimension of ambient spacetime
N_constr = 200; %Number of points on the circle at which constraints are applied
 
%alpha2 = (D-26)/(384*pi);
beta3 = 0;
 
dataOutFile=sprintf('Sm_data_D%d_n%d_N%d.dat',D,n,N_constr);
 
alpha2_values=[-0.1, (D-26)/384/pi, 0.1];%linspace(-0.2,0.2,100);
tv=zeros(1,length(alpha2_values));
count=1;
 
%fd = fopen(dataOutFile,'w');
 
for alpha2 = alpha2_values;
cvx_begin quiet;
    cvx_precision high;
    cvx_solver mosek;
    variable an(n);
    variable bn(n);
    variable cn(n);
     
     
    ns = (0:n-1)';
     
    %coefficients in s plane in terms of those in chi plane
    a0_t = sum(an);
    b0_t = sum(bn);
 
    a1_t = sum(-ns.*bn/2);
    b1_t =  sum(ns.*an/2);
 
    a2_t = sum(-ns.^2/8.*an);
    b2_t = sum(-ns.^2/8.*bn);
 
    a3_t = sum(ns.*(2*ns.^2+1)/96.*bn);
    b3_t = sum(-ns.*(2*ns.^2+1)/96.*an);
 
    c0_t = sum(cn);
 
    c1_t_by_i = sum(ns/2.*cn); %c1_t/1j
     
    c2_t = sum(-ns.^2/8.*cn);
 
    c3_t_by_i = sum(-ns.*(2*ns.^2+1)/96.*cn); %c3_t/1j
 
    t=c3_t_by_i-beta3;
 
    alpha3 = t+1/384
    minimize(t);
     
    a0_t == 0;
    b0_t == 0;
    a1_t == 0;
    b1_t == 0;
    a2_t == 0;
    b2_t == -alpha2;
    a3_t == alpha2/4;
    b3_t == -beta3;
    c0_t == 1;
    c1_t_by_i == 0.25;
    c2_t == -1/32;
     
    theta=linspace(0, pi, N_constr);
     
    ch = exp(1i*theta.'*(0:n-1));
    s_1 = ch*(an+1i*bn);
    s_2 = ch*cn;
    s_3 = ch*(an-1i*bn);
    s_sym= s_2+s_3;
    s_anti = s_2-s_3;
    s_sing = (D-2)*s_1+s_sym;
     
    real(s_sing).^2 + imag(s_sing).^2 <= 1;
    real(s_sym).^2  + imag(s_sym).^2  <= 1;
    real(s_anti).^2 + imag(s_anti).^2 <= 1;
      
      
 cvx_end
 
     
     fprintf('%4.20f %4.20f\n',alpha2,alpha3);
     %fprintf(fd,'%4.20f %4.20f\n',alpha2,alpha3);
     tv(count)=alpha3;
     count=count+1;
  
  
end
 
%fclose(fd);
%f = figure('visible','off');
%plot(alpha2_values,tv,'ro');
%saveas(f, 'bootstrap_beta3_zero.png','png')
