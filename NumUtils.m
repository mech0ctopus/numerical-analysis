classdef NumUtils
    % Defines numerical analysis utility functions.
    %
    % Methods:
    %   AB2 ------------- Adams-Bashforth 2-step (AB2) LMM for solving IVP.
    %   Bisection ------- Bisection Method.
    %   CentralDiff ----- Central Difference for approximating f'(x0).
    %   EstimateFPIter -- Estimates the number of iterations required for
    %                     convergence of Fixed Point Iteration algorithm.
    %   EulersMethod ---- Euler's Method for solving an IVP.
    %   FivePointMidpoint - Five-Point Midpoint Formula for approximating 
    %                       f'(x0).
    %   FixedPointIter -- Fixed Point Iteration.
    %   ForwardDiff ----- Forward Difference for approximating f'(x0).
    %   GaussSeidelMethod - Gauss-Seidel method for iteratively solving a
    %                       linear system of equations.
    %   GramSchmidt ----- Construct orthogonal polynomials w/ Gram-Schmidt.
    %   Jacobian -------- Symbolically calculates Jacobian for a system.
    %   JacobisMethod --- Jacobi's method for iteratively solving a
    %                     linear system of equations
    %   Lagrange -------- Generate Lagrange interpolating polynomial.
    %   LLS ------------- Constructs linear least squares poly. coeffs.
    %   LogB ------------ Calculates log(X) with base B.
    %   NaturalCubicSpline - Calculates the natural cubic spline for f.
    %   NewtonsMethod --- Newton's method for root finding problem
    %   NewtonsMethodForSystems - Newton's Method for iteratively solving a
    %                             nonlinear system of equations F(x)=0.
    %   QuasiNewton  ---- Quasi-Newton Method for root finding problem 
    %                     using global Bisection Method and local Newton's.
    %   QuasiSecant  ---- Quasi-Secant Method for root finding problem 
    %                     using global Bisection Method and local Secant.
    %   RK2 ------------- Runge-Kutta 2-step (RK2) for solving IVP.
    %   SecantMethod ---- Secant method for root finding problem
    %   TaylorPoly ------ Symbolically calculate the first N terms of a 
    %                     Taylor Polynomial.
    %   TaylorPolyNTerm - Symbolically calculate the Nth term of a Taylor
    %                     Polynomial.
    %   ThreePointEndpoint - Three-Point Endpoint Formula for approximating 
    %                        f'(x0).
    %   ThreePointMidpoint - Three-Point Midpoint Formula for approximating 
    %                        f'(x0).
    %   TruncationError - Symbolically calculate the truncation error 
    %                     associated with Taylor Polynomial approximation.
    %   TruncationErrorLagrange - Symbolically calculate the truncation 
    %                             error associated with Lagrange 
    %                             Interpolating Polynomial approximation.
    %
    % Usage: NumUtils.FunctionName(args)

    methods(Static)
        %% AB2
        function [t_vec,y_vec] = AB2(f,a,b,alpha,h,onestep)
            % Uses the Adams-Bashforth 2-step (AB2) linear multistep method
            % for solving an IVP.
            %
            % Input:  f = (y'(t) = f(t,y))
            %         a = Start point of time interval.
            %         b = End point of time interval.
            %         alpha = y(a)
            %         h = Timestep
            %         onestep = One-step method to use for computing y1.
            %                   Select forward_euler, backward_euler, or rk2.
            %
            % Output: t_vec = Column vector of mesh points (each timestep)
            %         y_vec = Approximated solution at each mesh point

            % Initialize output vectors
            t_vec=a:h:b; 
            y_vec=zeros(length(t_vec),1);
            y_vec(1)=alpha;

            % Use onestep method to get y_1
            switch onestep
                case 'forward_euler'
                    [~,w]=NumUtils.EulersMethod(f,a,a+h,h,alpha,'forward');
                    % Store result of y_1
                    y_vec(2)=w(2);
                    % Iteratively calculate solution at each mesh point
                    for ii=3:length(t_vec)
                        f_term=-f(t_vec(ii-2),y_vec(ii-2))+3*f(t_vec(ii-1),y_vec(ii-1));
                        y_vec(ii)=y_vec(ii-1)+(h./2).*f_term;
                    end
                case 'backward_euler'
                    [~,w]=NumUtils.EulersMethod(f,a,a+h,h,alpha,'backward');
                    % Store result of y_1
                    y_vec(2)=w(2);
                    % Iteratively calculate solution at each mesh point
                    for ii=3:length(t_vec)
                        f_term=-f(t_vec(ii-2),y_vec(ii-2))+3*f(t_vec(ii-1),y_vec(ii-1));
                        y_vec(ii)=y_vec(ii-1)+(h./2).*f_term;
                    end
                case 'rk2'
                    [w]= NumUtils.RK2(f,[a;a+h],alpha);
                    % Store result of y_1
                    y_vec(2)=w(2);
                    % Iteratively calculate solution at each mesh point
                    for ii=3:length(t_vec)
                        f_term=-f([t_vec(ii-2),y_vec(ii-2)])+3*f([t_vec(ii-1),y_vec(ii-1)]);
                        y_vec(ii)=y_vec(ii-1)+(h./2).*f_term;
                    end
                otherwise
                    disp('Error selecting onestep method. Please use ' ...
                         + 'forward_euler, backward_euler, or rk2.')
            end
        end
        %% Bisection
        function [p,iter,relerr,p_all,iter_all,relerr_all] = Bisection(f,a,b,tol,maxiter,verbose)
            % Function implementing the bisection method for solving the 
            % root-finding problem f(x) = 0 for a continuous function f on the 
            % closed interval [a,b]
            %
            % Need f(a) and f(b) to have different signs
            %
            % Input:  f   = @(x) function
            %         a   = left endpoint of x interval
            %         b   = right endpoint of x interval
            %         tol = tolerance for stopping criterion
            %         maxiter = maximum number of iterations
            %         verbose = prints extra information if equal to 1
            %
            % Output: p   = approximated root of f on [a,b]
            %         iter = total number of iterations performed
            %         relerr = resulting relative error
            %         p_all   = approximated root of f on [a,b] at each step
            %         iter_all = total number of iterations performed at each step
            %         relerr_all = resulting relative error at each step

            iter=0; FA=f(a);
            %Initialize output arrays
            p_all=[]; iter_all=[]; relerr_all=[];
            if verbose
                fprintf("Bisection Method Results:\n");
            end
            while iter<maxiter
                iter=iter+1;
                iter_all(iter)=iter;
                p=a+(b-a)/2; p_all(iter)=p;
                if verbose
                    fprintf('iter: %.f, a: %.2f, b: %.2f, p%.f: %.8f\n',iter,a,b,iter,p);
                end
                FP=f(p);
                relerr=(b-a)/2; relerr_all(iter)=relerr;
                % Use 1e-15 as tolerance when checking equal to zero since doubles 
                % are accurate to roughly 16 decimal places
                if ((FP<1e-15 && FP>=0) || (FP>-1e-15 && FP<=0) || relerr<tol)
                    break
                end


                if FA*FP>0
                    a=p; FA=FP;
                else
                    b=p;
                end
            end
            if verbose
                if iter>maxiter
                    fprintf('Method failed after %.0f iterations.\n', maxiter);
                end
            end
        end
        
        %% CentralDiff
        function [df_x0] = CentralDiff(f,h)   
            % Central-Difference Formula for approximating first-derivative 
            % at x0.
            %
            % Implemented per Equation 4.1 of Burden, Faires, Burden.
            %
            % Input:  f = Vector [f(x0+h) f(x0-h)]
            %         h = Distance between x-nodes.
            %
            % Output: df_x0 = Approximation of f'(x0)

            % Apply forward difference approximation
            df_x0 = 1/(2*h)*(f(1)-f(2));
        end
        
        %% EstimateFPIter
        function [N] = EstimateFPIter(g,k,p0,tol,verbose)
            % Function for estimating the number of iterations required to achieve
            % the desired tolerance using Fixed Point Iteration. Based on corollary
            % to Fixed-Point Theorem.
            %
            % Input:  f   = @(x) function
            %         k =  constant from Fixed-Point Theorem (0 < k < 1)
            %         p0 =  any number on the interval [a,b]
            %         tol = tolerance for stopping criterion
            %
            % Output: N   = estimated number of iterations to achieve tolerance

            % Calculate p1
            p1=g(p0);
            %Calculate left hand side of inequality and move terms over
            left_side=tol*(1-k)/abs(p1-p0);
            %Use log property: logB(X) = logA(X) / logA(B)
            N=log(left_side) / log(k);
            if verbose
                fprintf('Number of iterations required to converge ');
                fprintf('to fixed point within tolerance of %.8f is <=%.4f.\n',tol,N);
            end
        end
        
        %% EulersMethod
        function [tt,w]=EulersMethod(f,a,b,h,alpha,method)
            % Uses Euler's Method to approximate the solution of a well-posed IVP
            % at equally spaced numbers in the interval [a,b].
            %
            % Implemented per p.267 of Burden, Faires, Burden.
            %
            % Input:  f = anonymous function y'=f(t,y)
            %         a = start point of interval
            %         b = end point of interval
            %         h = step-size between mesh points
            %         alpha = initial condition y(a)=alpha
            %         method = difference method to use ('forward', 'midpoint')
            %
            % Output: tt = mesh points
            %         w  = approximation to y at each mesh point

            % Create time-vector (mesh points). 
            tt=a:h:b;
            % Initialize output vector
            w=zeros(length(tt),1); w(1)=alpha;

            switch method
                case 'forward'
                    % For each meshpoint, use forward approx.
                    for ii=2:length(tt)
                        % Calculate the y-approximation (w)
                        w(ii)=w(ii-1)+h.*f(tt(ii-1),w(ii-1));
                    end
                case 'backward'
                    % For each meshpoint, use backward approx.
                    tol=10^-6; maxiter=100; verbose=0;
                    for ii=2:length(tt)
                        % Calculate the y-approximation (w)
                        g=@(u) w(ii-1)+h.*f(tt(ii),u)-u;
                        % Use Bisection Method since backward-euler is implicit
                        [w(ii),~,~,~,~,~] = Bisection(g,w(ii-1)-1,w(ii-1)+1,...
                                                      tol,maxiter,verbose);
                    end
                case 'midpoint'
                    % Get second initial point using forward approx.
                    w(2)=w(1)+h.*f(tt(1),w(1));
                    % For each meshpoint afterwords, use midpoint method
                    for ii=3:length(tt)
                        % Calculate the y-approximation (w)
                        w(ii)=w(ii-2)+2.*h.*f(tt(ii-1),w(ii-1));
                    end
                otherwise
                    fprintf("Error, unknown method. ");
                    fprintf("Use 'forward' or 'midpoint'\n");
            end

        end
        %% FivePointMidpoint
        function [df_x0] = FivePointMidpoint(f,h)
            % Five-Point Endpoint Formula for approximating first-derivative 
            % at x0.
            %
            % Implemented per Equation 4.6 of Burden, Faires, Burden.
            %
            % Input:  f = Vector [f(x0-2*h) f(x0-h) f(x0+h) f(x0+2*h)]
            %         h = Distance between x-nodes.
            %
            % Output: df_x0 = Approximation of f'(x0)

            df_x0=1/(12*h).*(f(1)-8*f(2)+8*f(3)-f(4));
        end

        %% FixedPointIter
        function [p,iter,relerr] = FixedPointIter(g,p0,tol,maxiter)
            % Function implementing fixed-point iteration for g(x) on [a,b] to
            % find a solution to p=g(p) given an initial approximation p0.
            % Assumes g(x) has met existence and uniqueness criteria.
            %
            % Algorithm described in Numerical Analysis (Burden, Faires, Burden).
            %
            % Input:  f   = @(x) function
            %         p0 =  any number on the interval [a,b]
            %         tol = tolerance for stopping criterion
            %         maxiter = maximum number of iterations
            %
            % Output: p   = approximated unique fixed point of g on [a,b]
            %         iter = total number of iterations performed
            %         relerr = resulting relative error

            %Initialize values
            iter=1; relerr=inf;
            %Perform fixed point iteration
            while (iter<maxiter)
                %Compute p_i
                p=g(p0);
                %Check if error is within tolerance and procedure succeeded
                relerr=abs(p-p0);
                if relerr<tol
                    break
                end
                %Update parameters for next iteration
                iter=iter+1; p0=p;
            end

            %Display message if method failed
            if iter>=maxiter
                fprintf('Method failed after %d iterations.\n',iter);
            end
        end
        
        %% ForwardDiff
        function [df_x0] = ForwardDiff(f,h)   
            % Forward-Difference Formula for approximating first-derivative 
            % at x0.
            %
            % Implemented per Equation 4.1 of Burden, Faires, Burden.
            %
            % Input:  f = Vector [f(x0+h) f(x0)]
            %         h = Distance between x-nodes.
            %
            % Output: df_x0 = Approximation of f'(x0)

            % Apply forward difference approximation
            df_x0 = 1/h*(f(1)-f(2));
        end
        
        %% GaussSeidelMethod
        function [XO,iter,norm, ...
                  XO_all,iter_all,norm_all] = GaussSeidelMethod(A,b,x0,tol,maxiter,verbose)
            % Function implementing Gauss-Seidel method for iteratively solving a
            % linear system of equations.
            %
            % Implemented per p.461 of Burden, Faires, Burden.
            %
            % Input:  A   = A-matrix (from form Ax=b)
            %         b   = b-matrix (from form Ax=b)
            %         x0  = Initial approximation of x (from form Ax=b)
            %         tol = tolerance for stopping criterion
            %         maxiter = maximum number of iterations
            %         verbose = prints extra information if equal to 1
            %
            % Output: XO   = approximated solution for x (from form Ax=b)
            %         iter = total number of iterations performed
            %         norm = resulting infinity norm
            %         XO_all   = approximated solution for x at each step
            %         iter_all = total number of iterations performed at each step
            %         norm_all = resulting infinity norm at each step

            % Identify number of equations/unknowns
            n=length(b); 
            XO=x0; x=zeros(n,1);
            iter=0;
            %Initialize output arrays
            XO_all=[]; iter_all=[]; norm_all=[];
            if verbose
                fprintf("Gauss-Seidel Method Results:\n");
            end
            % Begin iteration
            while (iter<maxiter)
                % Update iter variable
                iter=iter+1;
                XO_all(iter,:)=XO; iter_all(iter,1)=iter;
                % Calculate new x-approximation
                for ii=1:n
                    % Calculate ax summation term for this step
                    ax_term=0;
                    for jj=1:(ii-1)
                        ax_term=ax_term+A(ii,jj).*x(jj,1);
                    end
                    % Calculate aXO summation term for this step
                    aXO_term=0;
                    for jj=(ii+1):n
                        aXO_term=aXO_term+A(ii,jj).*XO(jj,1);
                    end
                    x(ii,1)=(1./A(ii,ii)).*(-ax_term-aXO_term+b(ii,1));
                end
                % Check convergence with infinity norm
                norm=norm(x-XO,Inf); norm_all(iter,1)=norm;
                if verbose
                    fprintf('iter: %.f, norm: %.8f\n',iter,norm);
                end
                if norm<tol
                    break
                end
                % Update XO
                XO=x;
            end
            if verbose
                if iter>maxiter
                    fprintf('Method failed after %.0f iterations.\n', maxiter);
                end
            end
        end
        
        % GramSchmidt
        function [phi_k,Bk,Ck]=GramSchmidt(w,a,b,k,prev_phi,pprev_phi)
            % Uses Gram-Schmidt process to construct orthogonal polynomials
            % on the interval [a,b]
            %
            % Input:  w = weight function
            %         a = start point of interval
            %         b = end point of interval
            %         k = current k-value
            %         prev_phi  = k-1 polynomial function
            %         pprev_phi = k-2 polynomial function
            %
            % Output: phi_k = current (k) polynomial function (anonymous)
            %         Bk = B coefficient
            %         Ck = C coefficient

            % Calculate Bk
            B_num=@(x) x.*w.*(prev_phi(x).^2);
            B_den=@(x) w.*(prev_phi(x).^2)+0.*x;
            Bk=double(integral(B_num,a,b)./integral(B_den,a,b));

            if k==1
                Ck=0;
                % Define polynomial function for k==1
                phi_k=@(x) x-Bk;
            else
                % Calculate Ck for k>2
                C_num=@(x) x.*w.*prev_phi(x).*pprev_phi(x);
                C_den=@(x) w.*(pprev_phi(x).^2)+0.*x;
                Ck=double(integral(C_num,a,b)./integral(C_den,a,b));
                % Define polynomial function for k>2
                phi_k=@(x) (x-Bk).*prev_phi(x)-Ck.*pprev_phi(x);
            end
        end
        
        %% Jacobian
        function J=Jacobian(F,X)
            % Function for symbollically calculating the Jacobian for a system of
            % equations and unknown variables.
            %
            % Input:  F   = Column cell array of symbolic functions
            %         X   = Column vector of symbolic unknown variables
            %
            % Output: J   = Symbolic Jacobian matrix (n x n)

            % Get number of input functions
            n=length(F);
            % Make sure number of functions and unknowns match
            if n == length(X)
                % Initialize Jacobian as a symbolic matrix
                J=sym('j',[n n]);
                % For each row
                for ii=1:n
                    % For each column
                    for jj=1:n
                        % Calculate df(ii)/dx(jj)
                        J(ii,jj)=diff(F(ii),X(jj));
                    end
                end
            else
                disp('Number of functions and unknowns do not match.');
            end
        end
        
        %% JacobisMethod
        function [XO,iter,norm, ...
                  XO_all,iter_all,norm_all] = JacobisMethod(A,b,x0,tol,maxiter,verbose)
            % Function implementing Jacobi's method for iteratively solving a
            % linear system of equations.
            %
            % Implemented per p.459 of Burden, Faires, Burden.
            %
            % Input:  A   = A-matrix (from form Ax=b)
            %         b   = b-matrix (from form Ax=b)
            %         x0  = Initial approximation of x (from form Ax=b)
            %         tol = tolerance for stopping criterion
            %         maxiter = maximum number of iterations
            %         verbose = prints extra information if equal to 1
            %
            % Output: XO   = approximated solution for x (from form Ax=b)
            %         iter = total number of iterations performed
            %         norm = resulting infinity norm
            %         XO_all   = approximated solution for x at each step
            %         iter_all = total number of iterations performed at each step
            %         norm_all = resulting infinity norm at each step

            % Identify number of equations/unknowns
            n=length(b); 
            XO=x0; x=zeros(n,1);
            iter=0;
            %Initialize output arrays
            XO_all=[]; iter_all=[]; norm_all=[];
            if verbose
                fprintf("Jacobi's Method Results:\n");
            end
            % Begin iteration
            while (iter<maxiter)
                % Update iter variable
                iter=iter+1;
                XO_all(iter,:)=XO; iter_all(iter,1)=iter;
                % Calculate new x-approximation
                for ii=1:n
                    % Calculate summation term for this step
                    sum_term=0;
                    for jj=1:n
                        % Only include terms where jj is not equal to ii
                        if jj ~= ii
                            sum_term=sum_term+A(ii,jj).*XO(jj,1);
                        end
                    end
                    x(ii,1)=(1./A(ii,ii)).*(-sum_term+b(ii,1));
                end
                % Check convergence with infinity norm
                norm=norm(x-XO,Inf); norm_all(iter,1)=norm;
                if verbose
                    fprintf('iter: %.f, norm: %.8f\n',iter,norm);
                end
                if norm<tol
                    break
                end
                % Update XO
                XO=x;
            end
            if verbose
                if iter>maxiter
                    fprintf('Method failed after %.0f iterations.\n', maxiter);
                end
            end
        end
        
        %% Lagrange
        function poly = Lagrange(xpts,ypts,xeval)
            % Function to generate Lagrange interpolating polynomial at 
            % values xeval that passes through points (xpts,ypts)
            %
            % Input: xpts  = x points
            %        ypts  = y points
            %        xeval = evaluate polynomial at these x values
            %
            % Output: poly = Lagrange interpolating polynomial (order n)

            N = length(xpts); %n+1
            L = ones(N,length(xeval));
            poly = 0;

            % Generate Lagrange functions L_k(x) for each point xeval
            for k = 1:N
                for i = 1:N
                    if (i ~= k)
                        L(k,:) = L(k,:).*(xeval-xpts(i))./(xpts(k)-xpts(i));
                    end
                end
                % Generate Lagrange interpolating polynomial
                poly = poly + ypts(k)*L(k,:);
            end
        end
        
        %% LLS
        function [theta]=LLS(x,y,deg)
            % Function for constructing the linear least squares (LLS) polynomial
            % coefficients
            %
            % Input:  x   = Input datapoints/nodes, xn, column vector
            %         y   = Input datapoints, yn=f(xn), column vector
            %         deg = Degree of polynomial
            %
            % Output: theta = linear least squares polynomial coefficients

            % Construct design matrix
            X = zeros(length(x),deg);
            for ii=0:deg
                X(:,ii+1)=x.^ii;
            end
            % Solve for coefficients
            theta = (X'*X)\(X'*y);
            % Reverse order of theta for use with polyval (descending powers)
            theta=flip(theta);
        end   
        
        %% LogB
        function [logBX] = LogB(X,B)
            % Function for calculating log(x) with any base
            % using log property: logB(X) = logA(X) / logA(B)
            %
            % Input:  X   = value to take log of
            %         B =  Base
            %
            % Output: logBX   = log(X) with base 'B'

            logBX = log(X)/log(B);
        end
        
        %% NaturalCubicSpline
        function [a,b,c,d,S] = NaturalCubicSpline(x,y)
            % Function for constructing the cubic spline interpolant S for the
            % function f, defined at x0<x1<xn, satisfying S''(x0)=S''(xn)=0.
            %
            % Implemented per p.147 of Burden, Faires, Burden.
            %
            % Input:  x   = Input datapoints/nodes, xn, column vector
            %         y   = Input datapoints, yn=f(xn), column vector
            %
            % Output: a   = a-coefficients
            %         b   = b-coefficients
            %         c   = c-coefficients
            %         d   = d-coefficients
            %         S   = Symbolic cubic spline approximation of form:
            %               S(x)=sj(x)=aj+bj(x-xj)+cj(x-xj)^2+dj(x-xj)^3
            %               for xj<=x<=xj+1

            % Calculate number of data points. Define N.
            n=length(x)-1; N=n+1;
            % Define a-coefficients. Calculate h-values (step sizes between nodes)
            a=y; h = diff(x);

            % Build A matrix and b_vec (A*c=b_vec) to find c coefficients
            A=zeros(N,N);A(1,1)=1;A(end,end)=1; b_vec=zeros(N,1);
            for ii=2:n
                A(ii,ii-1)=h(ii-1);
                A(ii,ii)=2.*(h(ii-1)+h(ii));
                A(ii,ii+1)=h(ii);
                b_vec(ii,1)=3./h(ii).*(a(ii+1)-a(ii))-(3./h(ii-1)).*(a(ii)-a(ii-1));
            end

            % Solve for c
            c = A\b_vec;

            % Initialize b and d vectors;
            b=zeros(n,1);d=zeros(n,1);
            for jj=n:-1:1
                b(jj)=(a(jj+1)-a(jj))./h(jj)-h(jj).*(c(jj+1)+2.*c(jj))./3;
                d(jj)=(c(jj+1)-c(jj))./(3.*h(jj));
            end         
            a=a(1:n,1); c=c(1:n,1);

            % Calculate S symbollically
            S=sym('x',[n 1]); syms xx real;
            for jj=n:-1:1
                S(jj)=a(jj)+b(jj).*(xx-x(jj))+c(jj).*(xx-x(jj)).^2+d(jj).*(xx-x(jj)).^3;
            end
        end

        %% NewtonsMethod
        function [p,iter,relerr, ...
                  p_all,iter_all,relerr_all] = NewtonsMethod(f,p0,tol,maxiter, ...
                                                             verbose)
            % Function implementing Newton's method for solving the root-finding
            % problem f(x) = 0 given an initial approximation p0
            %
            % Implemented per p.67 of Numerical Analysis by Burden, Faires, Burden.
            %
            % Input:  f   = @(x) function
            %         p0 = Initial approximation to p in [a,b]
            %         tol = tolerance for stopping criterion
            %         maxiter = maximum number of iterations
            %         verbose = prints extra information if equal to 1
            %
            % Output: p   = approximated root of f
            %         iter = total number of iterations performed
            %         relerr = resulting relative error
            %         p_all   = approximated root of f on [a,b] at each step
            %         iter_all = total number of iterations performed at each step
            %         relerr_all = resulting relative error at each step

            % Calculate derivative of f
            syms x; df=diff(f,x);
            % Perform Newton's Method
            if verbose
                fprintf("Newton's Method Results:\n");
                fprintf('iter: 0, p0: %.16f,\n',p0);
            end
            iter=0;
            %Initialize output arrays
            p_all=[]; iter_all=[]; relerr_all=[];
            while iter<maxiter
                % Update iteration counter 
                iter=iter+1; iter_all(iter)=iter;
                % Evaluate df(p0) and get numerical result
                df_p0=double(subs(df,p0));
                % Calculate p
                p=p0-f(p0)./df_p0; p_all(iter)=p;
                % Print information about current iteration
                if verbose
                    fprintf('iter: %.f, p%.f: %.16f\n', iter, iter, p);
                end
                % Check convergence
                relerr=abs(p-p0); relerr_all(iter)=relerr;
                if relerr<tol
                    break
                end
                % Update p0
                p0=p;
            end

            if iter>maxiter
                fprintf('Method failed after %.0f iterations', maxiter);
            end
        end   
        
        %% NewtonsMethodForSystems
        function [XO,iter,norm, ...
                  XO_all,iter_all,norm_all] = NewtonsMethodForSystems(F,X,x0,tol,maxiter,verbose)
            % Function implementing Newton's Method for iteratively solving a
            % nonlinear system of equations F(x)=0.
            %
            % Implemented per p.653 of Burden, Faires, Burden.
            %
            % Input:  F   = Column cell array of nonlinear mappings
            %         X   = Column vector of symbolic variables
            %         x0  = Initial approximation of x
            %         tol = tolerance for stopping criterion
            %         maxiter = maximum number of iterations
            %         verbose = prints extra information if equal to 1
            %
            % Output: XO   = approximated solution for x (from form Ax=b)
            %         iter = total number of iterations performed
            %         norm = resulting infinity norm
            %         XO_all   = approximated solution for x at each step
            %         iter_all = total number of iterations performed at each step
            %         norm_all = resulting infinity norm at each step

            % Initalize variables
            XO=x0; iter=0;
            % Initialize output arrays
            XO_all=[]; iter_all=[]; norm_all=[];
            if verbose
                fprintf("Newton's Method for Systems Results:\n");
            end
            % Begin iteration
            while (iter<maxiter)
                % Update iter variable
                iter=iter+1;
                XO_all(iter,:)=XO; iter_all(iter,1)=iter;
                % Calculate F_x by substituting x-values in and converting to
                % double
                F_x=subs(F,'x1',XO(1));
                F_x=subs(F_x,'x2',XO(2));
                F_x=double(subs(F_x,'x3',XO(3)));

                % Calculate symbolic Jacobian
                J=NumUtils.Jacobian(F,X);
                % Calculate J_x by substituting x-values in and converting to
                % double
                J_x=subs(J,'x1',XO(1));
                J_x=subs(J_x,'x2',XO(2));
                J_x=double(subs(J_x,'x3',XO(3)));

                % Solve linear system (n x n)
                % If the Jacobian is singular/non-invertible
                if cond(J_x)==Inf
                    disp('J(x) is singular.');
                    break
                else
                    y=J_x\-F_x;
                end
                % Update x
                XO=XO+y;
                % Check convergence with infinity norm
                norm=norm(y,Inf); norm_all(iter,1)=norm;
                if verbose
                    fprintf('iter: %.f, norm: %.8f\n',iter,norm);
                end
                if norm<tol
                    break
                end
            end
            if verbose
                if iter>maxiter
                    fprintf('Method failed after %.0f iterations.\n', maxiter);
                end
            end
        end
        
        %% QuasiNewton
        function [p,iter,relerr, ...
                  p_all,iter_all, ...
                  relerr_all]=QuasiNewton(f,a,b,tol,max_global_steps, ...
                                                     max_local_steps,verbose)
            % Function implementing Quasi-Newton scheme for solving the root-finding
            % problem f(x) = 0 for a continuous function f on the closed interval 
            % [a,b]. Uses Bisection as global method to approximate p0, then
            % Newton's Method for final convergence.
            %
            % Input:  f   = @(x) function
            %         a   = left endpoint of x interval
            %         b   = right endpoint of x interval
            %         tol = tolerance for stopping criterion
            %         max_global_steps = maximum number of Bisection iterations
            %         max_local_steps = maximum number of Newton's Method iterations
            %         verbose = prints extra information if equal to 1
            %
            % Output: p   = approximated root of f on [a,b]
            %         iter = total number of iterations performed
            %         relerr = resulting relative error
            %         p_all   = approximated root of f on [a,b] at each step
            %         iter_all = total number of iterations performed at each step
            %         relerr_all = resulting relative error at each step
            %         iter_type = identify whether global or local step

            % Perform global bisection steps to get a good p0 approximation
            [p0,iter_global,relerr_global, ...
             p0_all,iter_global_all, ...
             relerr_global_all] = Bisection(f,a,b,tol,max_global_steps,0);

            % Store global results in final output arrays
            p_all=p0_all; iter_all=iter_global_all; relerr_all=relerr_global_all;
            iter_type(1:length(p0_all))="global";

            % Perform local Newton's Method steps until convergence
            [p,iter,relerr, ...
             p_local_all,iter_local_all, ...
             relerr_local_all] = NewtonsMethod(f,p0,tol,max_local_steps,0);

            % Add local results to final output arrays
            iter_local_all=iter_local_all+iter_global;
            p_all=[p_all,p_local_all]; iter_all=[iter_all,iter_local_all]; 
            relerr_all=[relerr_all,relerr_local_all];
            iter_type(length(p0_all)+1:length(p_all))="local";

            if verbose
                fprintf("Quasi-Newton Results:\n");
                for ii=1:length(p_all)
                    fprintf('iter: %.f, ',iter_all(ii));
                    fprintf('p%.f: %.8f, ',iter_all(ii),p_all(ii));
                    fprintf('relerr: %.4f, ', relerr_all(ii));
                    fprintf('type: %s\n', iter_type(ii));
                end
            end    

        end
        
        %% QuasiSecant
        function [p,iter,relerr, ...
                  p_all,iter_all, ...
                  relerr_all]=QuasiSecant(f,a,b,tol,max_global_steps, ...
                                                     max_local_steps,verbose)
            % Function implementing Quasi-Secant scheme for solving the root-finding
            % problem f(x) = 0 for a continuous function f on the closed interval 
            % [a,b]. Uses Bisection as global method to approximate p0 and p1, then
            % Secant Method for final convergence.
            %
            % Input:  f   = @(x) function
            %         a   = left endpoint of x interval
            %         b   = right endpoint of x interval
            %         tol = tolerance for stopping criterion
            %         max_global_steps = maximum number of Bisection iterations
            %         max_local_steps = maximum number of Secant Method iterations
            %         verbose = prints extra information if equal to 1
            %
            % Output: p   = approximated root of f on [a,b]
            %         iter = total number of iterations performed
            %         relerr = resulting relative error
            %         p_all   = approximated root of f on [a,b] at each step
            %         iter_all = total number of iterations performed at each step
            %         relerr_all = resulting relative error at each step
            %         iter_type = identify whether global or local step

            % Perform global bisection steps to get a good p0 and p1 approximation
            [p_global,iter_global,relerr_global, ...
             p_global_all,iter_global_all, ...
             relerr_global_all] = Bisection(f,a,b,tol,max_global_steps,0);
            % Index results for p0 and p1
            p0=p_global_all(length(p_global_all)-1);
            p1=p_global_all(length(p_global_all));

            % Store global results in final output arrays
            p_all=p_global_all; iter_all=iter_global_all; 
            relerr_all=relerr_global_all;
            iter_type(1:length(p_global_all))="global";

            % Perform local Secant Method steps until convergence
            [p,iter,relerr, ...
             p_local_all,iter_local_all, ...
             relerr_local_all] = SecantMethod(f,p0,p1,tol, ...
                                              max_local_steps,0);

            % Add local results to final output arrays
            iter_local_all=iter_local_all+iter_global;
            p_all=[p_all,p_local_all]; iter_all=[iter_all,iter_local_all]; 
            relerr_all=[relerr_all,relerr_local_all];
            iter_type(length(p_global_all)+1:length(p_all))="local";

            if verbose
                fprintf("Quasi-Secant Results:\n");
                for ii=1:length(p_all)
                    fprintf('iter: %.f, ',iter_all(ii));
                    fprintf('p%.f: %.8f, ',iter_all(ii),p_all(ii));
                    fprintf('relerr: %.4f, ', relerr_all(ii));
                    fprintf('type: %s\n', iter_type(ii));
                end
            end    

        end
        
        %% SecantMethod
        function [p,iter,relerr, ...
                  p_all,iter_all,relerr_all] = SecantMethod(f,p0,p1,tol, ...
                                                            maxiter,verbose)
            % Function implementing the Secant method for solving the root-finding
            % problem f(x) = 0 given initial approximations p0 and p1.
            %
            % Implemented per p.71 of Numerical Analysis by Burden, Faires, Burden.
            %
            % Input:  f   = @(x) function
            %         p0 = Initial approximation to p in [a,b]
            %         p1 = Initial approximation of p1
            %         tol = tolerance for stopping criterion
            %         maxiter = maximum number of iterations
            %         verbose = prints extra information if equal to 1
            %
            % Output: p   = approximated root of f
            %         iter = total number of iterations performed
            %         relerr = resulting relative error
            %         p_all   = approximated root of f on [a,b] at each step
            %         iter_all = total number of iterations performed at each step
            %         relerr_all = resulting relative error at each step

            % Perform Secant Method
            iter=1; q0=f(p0); q1=f(p1);
            if verbose
                fprintf("Secant Method Results:\n");
                fprintf('iter: 0, p0: %.16f\n', p0);
                fprintf('iter: 1, p1: %.16f\n', p1);
            end
            %Initialize output arrays
            p_all=[p1]; iter_all=[1]; relerr_all=[abs(p1-p0)];
            while iter<maxiter
                % Update iteration counter
                iter=iter+1; iter_all(iter)=iter;
                % Calculate p
                p=p1-q1.*(p1-p0)./(q1-q0); p_all(iter)=p;
                % Print information about current iteration
                if verbose
                    fprintf('iter: %.f, p%.f: %.16f\n', iter, iter, p)
                end
                % Check convergence
                relerr=abs(p-p1); relerr_all(iter)=relerr;
                if relerr<tol
                    break
                end
                %Update p0, q0, p1, q1
                p0=p1; q0=q1; p1=p; q1=f(p);
            end

            if iter>maxiter
                fprintf('Method failed after %.0f iterations', maxiter);
            end
        end
        
        %% RK2
        function y = RK2(f,t,y0)
            % Function to compute RK2 approximation of solution to IVP y'=f(t,y) with
            % initial value y0 for t over the interval [a,b]
            %
            % Input: f  = RHS function of DE
            %        t  = mesh points of t values over the interval [a,b]
            %        y0 = initial value
            %
            % Output: y = RK2 approximation of solution to IVP

            h = t(end)-t(end-1);  % step size
            N = length(t);  % number of points in mesh / approximation

            y = NaN(N,1);
            y(1) = y0;

            for i = 1:N-1
                t_half = t(i)+h/2; % time points halway between mesh points
                y(i+1) = y(i) + h*f([t_half,y(i)+h/2*f([t(i),y(i)])]);  % RK2 formula
            end
        end
        
        %% TayloyPoly
        function [Pn]=TaylorPoly(f,x0,N)
            % Function for symbolically calculating the first N-terms of a
            % Taylor Polynomial.
            %
            % Input:  f   = anonymous @(x) function to approximate
            %         x0   = centerpoint of Taylor series
            %         N   = number of terms to calculate
            %
            % Output: Pn   = symbolic first N-terms of Taylor Polynomial 
            %                w.r.t. x.

            syms x;
            Pn=0;
            for k=0:N
                %Get kth derivative of f
                fk=diff(f,x,k);
                %Calculate current term of polynomial
                Pn_k=subs(fk,x0)/factorial(k)*((x-x0)^k);
                %Add current term to total
                Pn=Pn+Pn_k;
            end
        end
        
        %% TaylorPolyNTerm
        function [PN]=TaylorPolyNTerm(f,x0,N)   
            % Function for symbolically calculating the Nth term only of a
            % Taylor Polynomial.
            %
            % Input:  f   = anonymous @(x) function to approximate
            %         x0   = centerpoint of Taylor series
            %         N   = number of term to calculate
            %
            % Output: PN   = symbolic Nth term of Taylor Polynomial w.r.t. x.

            syms x;
            %Get Nth derivative of f
            fN=diff(f,x,N);
            %Calculate current term of polynomial
            PN=subs(fN,x0)/factorial(N)*((x-x0)^N);
        end
        
        %% ThreePointEndpoint
        function [df_x0] = ThreePointEndpoint(f,h)
            % Three-Point Endpoint Formula for approximating first-derivative 
            % at x0.
            %
            % Implemented per Equation 4.4 of Burden, Faires, Burden.
            %
            % Input:  f = Vector [f(x0) f(x0+h) f(x0+2*h)]
            %         h = Distance between x-nodes.
            %
            % Output: df_x0 = Approximation of f'(x0)

            df_x0=1/(2*h).*(-3*f(1)+4*f(2)-f(3));
        end
        
        %% ThreePointMidpoint
        function [df_x0] = ThreePointMidpoint(f,h)
            % Three-Point Midpoint Formula for approximating first-derivative 
            % at x0.
            %
            % Implemented per Equation 4.5 of Burden, Faires, Burden.
            %
            % Input:  f = Vector [f(x0+h) f(x0-h)]
            %         h = Distance between x-nodes.
            %
            % Output: df_x0 = Approximation of f'(x0)

            df_x0=1/(2*h).*(f(1)-f(2));
        end
        
        %% TruncationError
        function [Rn]=TruncationError(f,x0,n,ksi_x)
            % Function for symbolically calculating the truncation error 
            % associated with Taylor Polynomial approximation.
            %
            % Input:  f   = anonymous @(x) function to approximate
            %         x0   = centerpoint of Taylor series
            %         n   = order of Taylor Polynomial
            %         ksi_x = unknown function between x0 and x
            %
            % Output: Rn   = symbolic Truncation Error w.r.t. x.

            syms x;
            %Get (n+1)th derivative of f
            f_n1=diff(f,x,n+1);
            %Calculate truncation error
            Rn=subs(f_n1,ksi_x)/factorial(n+1)*((x-x0)^(n+1));
        end
        
        %% TruncationErrorLagrange
        function [Rn]=TruncationErrorLagrange(f,xn,ksi_x)
            % Function for symbolically calculating the truncation error 
            % associated with Lagrange Interpolating Polynomial approximation.
            %
            % Input:  f   = anonymous @(x) function to approximate
            %         xn   = vector of x-points used for interpolation
            %         ksi_x = value of unknown function (between x0 and xn)
            %
            % Output: Rn   = symbolic Truncation Error w.r.t. x.

            syms x;
            %Get (n+1)th derivative of f
            n=length(xn);
            f_n1=diff(f,x,n+1);
            %Calculate truncation error
            Rn=subs(f_n1,ksi_x)/factorial(n+1);
            for ii=1:n
                Rn=Rn*(x-xn(ii));
            end
        end
    end
end