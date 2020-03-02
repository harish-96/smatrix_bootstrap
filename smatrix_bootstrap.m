n = 50; %Order of taylor expansion
D = 4; %Dimension of ambient spacetime
N_constr = 200; %Number of points on the circle at which constraints are applied

alpha2 = (D-26)/(384*pi);
beta3 = 0;
cvx_begin sdp
    variable an(n);
    variable bn(n);
    variable cn(n);
    
    expression s_sym(N_constr);
    expression s_anti(N_constr);
    expression s_sing(N_constr);
    expression s_1(N_constr);
    expression s_2(N_constr);
    expression s_3(N_constr);
    expression M_sing(2,2,N_constr);
    expression M_sym(2,2,N_constr);
    expression M_anti(2,2,N_constr);
    
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

    minimize(c3_t_by_i);
    
    subject to
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
    
    i = 1;
     for theta=linspace(0, pi, N_constr)
        chi = exp(1j*theta);
        ch = chi.^(0:n-1);
        s_1(i) = sum(dot(an+1j*bn,ch));
        s_2(i) = sum(dot(cn,ch));
        s_3(i) = sum(dot(an-1j*bn,ch));
        s_sym(i) = s_2(i) + s_3(i);
        s_anti(i) = s_2(i) - s_3(i);
        s_sing(i) = (D-2)*s_1(i) + s_sym(i);
        
        M_sing(1,1,i) = 1+real(s_sing(i));
        M_sing(2,2,i) = 1-real(s_sing(i));
        M_sing(1,2,i) = imag(s_sing(i));
        M_sing(2,1,i) = imag(s_sing(i));
        
        M_sym(1,1,i) = 1+real(s_sym(i));
        M_sym(2,2,i) = 1-real(s_sym(i));
        M_sym(1,2,i) = imag(s_sym(i));
        M_sym(2,1,i) = imag(s_sym(i));
        
        M_anti(1,1,i) = 1+real(s_anti(i));
        M_anti(2,2,i) = 1-real(s_anti(i));
        M_anti(1,2,i) = imag(s_anti(i));
        M_anti(2,1,i) = imag(s_anti(i));
        
        M_sing(:,:,i) >= 0;
        M_anti(:,:,i) >= 0;
        M_sym(:,:,i) >= 0;
        
%         abs(s_sing(i)) <= 1;
%         abs(s_sym(i)) <= 1;
%         abs(s_anti(i)) <= 1;
%         real(s_sing(i))^2 + imag(s_sing(i))^2 <= 1;
%         real(s_sym(i))^2 + imag(s_sym(i))^2 <= 1;
%         real(s_anti(i))^2 + imag(s_anti(i))^2 <= 1;
        i = i+1;
     end

     
 cvx_end
