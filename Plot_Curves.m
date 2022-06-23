clear;
close all;
clc;
%% Get type of curve
type = questdlg('Which curve do you want?', 'Select curve type', 'Hermite','Bezier', 'Spline', 'default-nothing');

%% Get parameters & Plot
switch type
    
    case 'Hermite'
        u=0:0.01:1;
        parameter_titles = {'Enter p1_x:','Enter p1_y:','Enter p1''_x:','Enter p1''_y:','Enter p2_x:','Enter p2_y:','Enter p2''_x:','Enter p2''_y:'};
        parameters = num2cell(str2double(inputdlg(parameter_titles,'Enter the parameters',1)));
        [p_x, p_y, pr_x, pr_y, q_x, q_y, qr_x, qr_y]= parameters{:};
        x = (1-3*u.^2+2*u.^3)*p_x + (3*u.^2-2*u.^3)*q_x + (u-2*u.^2+u.^3)*pr_x + (-u.^2+u.^3)*qr_x;
        y = (1-3*u.^2+2*u.^3)*p_y + (3*u.^2-2*u.^3)*q_y + (u-2*u.^2+u.^3)*pr_y + (-u.^2+u.^3)*qr_y;
        
        plot(x,y)
        grid on
        hold on
        scatter([p_x; q_x], [p_y; q_y],'fill')
        
    case 'Bezier'
        order = questdlg('3D or 2D?', ...
        'Select curve order', '3D','2D','default_nothing');
    
    if order=='2D'
        n_str = inputdlg('Enter the degree (you need n+1 control points for n degree Bezier curve):','Degree of Bezier Curve',1);
        n = str2double(n_str);
        p = zeros(n+1,2);
        
        for i=1:(n+1)
            control_points_title = {['x-coordinate of point #',num2str(i)],['y-coordinate of point #',num2str(i)]};
            control_points = num2cell(str2double(inputdlg(control_points_title,'Enter the control points'' coordinates',1)));
            [p(i,1), p(i,2)] = control_points{:};
        end
        
        syms t;
        b = bernsteinMatrix(n,t);
        BezierCurve = simplify(b*p);
        fplot(BezierCurve(1),BezierCurve(2),[0,1])
        hold on
        grid on
        scatter(p(:,1), p(:,2),'fill')
        xlabel('x')
        ylabel('y')
        title(['Bezier curve of degree', ' ', num2str(n)])
        
    elseif order=='3D'
        n_str = inputdlg('Enter the degree (you need n+1 control points for n degree Bezier curve):','Degree of Bezier Curve',1);
        n = str2double(n_str);
        p = zeros(n+1,3);
        
        for i=1:(n+1)
            control_points_title = {['x-coordinate of point #',num2str(i)],['y-coordinate of point #',num2str(i)],['z-coordinate of point #',num2str(i)]};
            control_points = num2cell(str2double(inputdlg(control_points_title,'Enter the control points'' coordinates',1)));
            [p(i,1), p(i,2), p(i,3)] = control_points{:};
        end
        
        syms t;
        b = bernsteinMatrix(n,t);
        BezierCurve = simplify(b*p);
        fplot3(BezierCurve(1),BezierCurve(2),BezierCurve(3),[0,1])
        hold on
        grid on
        scatter3(p(:,1), p(:,2), p(:,3),'fill')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title(['Bezier curve of degree', ' ', num2str(n)])
    end
    
    case 'Spline'
    
    splinetype = questdlg('Which kind of spline do you want?', 'Select spline type', 'B-Spline', 'Cubic Curve', 'default-nothing');
    
    switch splinetype
        
    case 'B-Spline'
        
    order = questdlg('3D or 2D?', ...
    'Select curve order', '3D', '2D', 'default_nothing');
    

    if order=='2D'
        dp_str = inputdlg('Enter the degree of the curve (from zero to infinity):','Degree of B-Spline Curve',1);
        n_str = inputdlg('Number of control points - 1 (n):','Enter the number of control points "minus one"',1);
        dp = str2double(dp_str);
        n = str2double(n_str);
        p = zeros(n+1,2);
        
        for i=1:(n+1)
            control_points_title = {['x-coordinate of point #',num2str(i)],['y-coordinate of point #',num2str(i)]};
            control_points = num2cell(str2double(inputdlg(control_points_title,'Enter the control points'' coordinates',1)));
            [p(i,1), p(i,2)] = control_points{:};
        end
        
        
        % Definition of Basis Function in bspline_basis.m
        
        t = linspace(0,1,n+dp+1); %for uniformly distributed knots
        BSpline_curve_x = 0;
        BSpline_curve_y = 0;
        for i=1:(n+1)
            BSpline_curve_x = BSpline_curve_x + bspline_basis(i-1,dp,t)*p(i,1);
            BSpline_curve_y = BSpline_curve_y + bspline_basis(i-1,dp,t)*p(i,2);
        end
                        
        plot(BSpline_curve_x,BSpline_curve_y)
        hold on
        grid on
        scatter(p(:,1), p(:,2),'fill')
        xlabel('x')
        ylabel('y')
        title(['B-Spline curve of degree ', num2str(dp), ' ', 'with ', num2str(n), ' ', 'control points'])
        
    elseif order=='3D'
        dp_str = inputdlg('Enter the degree of the curve (from zero to infinity):','Degree of B-Spline Curve',1);
        n_str = inputdlg('Number of control points - 1 (n):','Enter the number of control points "minus one"',1);
        dp = str2double(dp_str);
        n = str2double(n_str);
        p = zeros(n+1,3);
        
        for i=1:(n+1)
            control_points_title = {['x-coordinate of point #',num2str(i)],['y-coordinate of point #',num2str(i)],['z-coordinate of point #',num2str(i)]};
            control_points = num2cell(str2double(inputdlg(control_points_title,'Enter the control points'' coordinates',1)));
            [p(i,1), p(i,2), p(i,3)] = control_points{:};
        end
        
        
        % Definition of Basis Function in bspline_basis.m
        
        t = linspace(0,1,n+dp+1); %for uniformly distributed knots
        BSpline_curve_x = 0;
        BSpline_curve_y = 0;
        BSpline_curve_z = 0;
        for i=1:(n+1)
            BSpline_curve_x = BSpline_curve_x + bspline_basis(i-1,dp,t)*p(i,1);
            BSpline_curve_y = BSpline_curve_y + bspline_basis(i-1,dp,t)*p(i,2);
            BSpline_curve_z = BSpline_curve_z + bspline_basis(i-1,dp,t)*p(i,3);
        end
                        
        plot3(BSpline_curve_x,BSpline_curve_y,BSpline_curve_z)
        hold on
        grid on
        scatter3(p(:,1), p(:,2), p(:,3),'fill')
        xlabel('x')
        ylabel('y')
        ylabel('z')
        title(['3D B-Spline curve of degree ', num2str(dp), ' ', 'with ', num2str(n), ' ', 'control points'])
    end
    
        case 'Cubic Curve'

        order = questdlg('3D or 2D?', ...
        'Select curve order', '3D','2D','default_nothing');
    
        n_str = inputdlg('Enter the number of control points','Number of Control Points',1);
        n = str2double(n_str);
        p = zeros(n,3);
                
        if order == '2D'
        for i=1:n
            control_points_title = {['x-coordinate of point #',num2str(i)],['y-coordinate of point #',num2str(i)]};
            control_points = num2cell(str2double(inputdlg(control_points_title,'Enter the control points'' coordinates',1)));
            [p(i,1), p(i,2)] = control_points{:};
        end
        
        syms u
        b = [(1-5.5*u+9*u^2-4.5*u^3), (9*u-22.5*u^2+13.5*u^3), (-4.5*u+18*u^2-13.5*u^3), (u-4.5*u^2+4.5*u^3)];
        
        CubicCurve = simplify(b*p(1:4,:));

        j=2;
        while (i~=0 && n~=4)
            if (j*4-1)<n
                CubicCurve = CubicCurve + simplify(b*p((j-1)*4:j*4-1,:));
            elseif (j*4-1)>n
                CubicCurve = CubicCurve + simplify(b(1:n-(j-1)*4+1)*p(4+(j-2)*3:n,:));
                i=0;
            end
            j=j+1;
        end
        fplot(CubicCurve(1),CubicCurve(2),[0,1])
        hold on
        grid on
        scatter(p(:,1), p(:,2),'fill')
        xlabel('x')
        ylabel('y')
        title(['Cubic Spline of degree', ' ', num2str(n)])
        
        elseif order == '3D'
        for i=1:(n)
            control_points_title = {['x-coordinate of point #',num2str(i)],['y-coordinate of point #',num2str(i)],['z-coordinate of point #',num2str(i)]};
            control_points = num2cell(str2double(inputdlg(control_points_title,'Enter the control points'' coordinates',1)));
            [p(i,1), p(i,2), p(i,3)] = control_points{:};
        end
        
        syms u
        b = [(1-5.5*u+9*u^2-4.5*u^3), (9*u-22.5*u^2+13.5*u^3), (-4.5*u+18*u^2-13.5*u^3), (u-4.5*u^2+4.5*u^3)];
        
        CubicCurve = simplify(b*p(1:4,:));

        j=2;
        while (i~=0 && n~=4)
            if (j*4-1)<n
                CubicCurve = CubicCurve + simplify(b*p((j-1)*4:j*4-1,:));
            elseif (j*4-1)>n
                CubicCurve = CubicCurve + simplify(b(1:n-(j-1)*4+1)*p(4+(j-2)*3:n,:));
                i=0;
            end
            j=j+1;
        end
        fplot3(CubicCurve(1),CubicCurve(2),CubicCurve(3),[0,1])
        hold on
        grid on
        scatter3(p(:,1), p(:,2),p(:,3),'fill')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title(['3D Cubic Spline curve of degree', ' ', num2str(n)])
        
        end
        
    end
        
end