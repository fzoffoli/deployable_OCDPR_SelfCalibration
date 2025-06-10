classdef RecordType
    properties
        n_cables;
        figure_handle;
        axes_handle;
        lines;
        frame_vector;
    end
    methods
        function obj = RecordType(par,title,name,CableForceCnt)
            obj.n_cables = par.n_cables;
            obj = obj.SetFigureParameter(par,title,name);
            for i=1:obj.n_cables
                obj.lines.pulley(i) = animatedline('Color','r','LineWidth',2);
                if i==CableForceCnt(1)||i==CableForceCnt(end) %%% FZ added this condition to show the force controlled cables
                    obj.lines.cable(i) = animatedline('Color','#D95319','LineWidth',2);
                else
                    obj.lines.cable(i) = animatedline('Color','k','LineWidth',2);
                end
                obj.lines.platformPoint(i) = animatedline('Color','b','LineWidth',1);
                obj.lines.platformEdge(i) = animatedline('Color','b','LineWidth',1);
                obj.lines.platformPlane(i) = animatedline('Color','b','LineWidth',1);
                obj.lines.PulleyAxis(i) = animatedline('Color','#D95319','LineStyle','-.','LineWidth',2);
            end
        end
        function obj = SetFigureParameter(obj,par,title,name)
            obj.figure_handle = figure('Name',title,'FileName',name,'NumberTitle','off','Position',[0 0 700 700]);
            for i=1:par.n_cables
                if (norm(par.cable(i,1).pos_OA_glob)>norm(par.cable(i,1).pos_PD_loc))
                    limits(:,i) = par.cable(i,1).pos_OA_glob;
                else
                    limits(:,i) = par.cable(i,1).pos_PD_loc;
                end
            end
            hold on
            max_x = max(limits(1,:)); min_x = min(limits(1,:));
            if (max_x==min_x)
                max_x = max_x+0.01;
                min_x = min_x-0.01;
            end
            incr_x = 0.02*(max_x - min_x);
            lim_x = [((min_x-incr_x)*10) ((max_x+incr_x)*10)]./10; 
            l_x = (lim_x(2)-lim_x(1)); tickSpace_x = l_x/6;
            max_y = max(limits(2,:)); min_y = min(limits(2,:));
            if (max_y==min_y)
                max_y = max_y+0.01;
                min_y = min_y-0.01;
            end
            incr_y = 0.02*(max_y - min_y);
            lim_y = [((min_y-incr_y)*10) ((max_y+incr_y)*10)]./10; l_y = (lim_y(2)-lim_y(1)); tickSpace_y = l_y/2;
            max_z = max(limits(3,:)); min_z = min(limits(3,:)); 
            if (max_z==min_z)
                max_z = max_z+0.01;
                min_z = min_z-0.01;
            end
            incr_z = 0.02*(max_z - min_z);
            middle_z = (max_z - min_z)/2;
            max_l = max([lim_x lim_y]); lim_z =[((min_z-incr_z)*10) ((max_z+incr_z)*10)]./10; l_z = (lim_z(2)-lim_z(1));
            tickSpace_z = l_z/6;
            
            obj.axes_handle = gca; obj.axes_handle.Box = 'Off'; obj.axes_handle.LineWidth = 2;
            obj.axes_handle.XAxisLocation = 'origin'; obj.axes_handle.XLim = lim_x;
            obj.axes_handle.XTick = lim_x(1):tickSpace_x:lim_x(2); obj.axes_handle.XGrid = 'On';
            obj.axes_handle.YAxisLocation = 'origin'; obj.axes_handle.ZLim = lim_z;
            obj.axes_handle.ZTick = lim_z(1):tickSpace_z:lim_z(2);   obj.axes_handle.YGrid = 'On';
            obj.axes_handle.YLim = lim_y; obj.axes_handle.YTick = lim_y(1):tickSpace_y:lim_y(2);
            obj.axes_handle.ZGrid = 'On'; obj.axes_handle.GridLineStyle = '--';%'-'
            
            obj.axes_handle.XLabel.String = 'x [m]'; obj.axes_handle.XLabel.FontSize = 16;
            obj.axes_handle.XLabel.FontWeight = 'bold';
            obj.axes_handle.YLabel.String = 'y [m]'; obj.axes_handle.YLabel.FontSize = 16;
            obj.axes_handle.YLabel.Rotation = 0; obj.axes_handle.YLabel.FontWeight = 'bold';
            obj.axes_handle.ZLabel.String = 'z [m]'; obj.axes_handle.ZLabel.FontSize = 16;
            obj.axes_handle.ZLabel.Rotation = 0; obj.axes_handle.ZLabel.FontWeight = 'bold';
            
            obj.axes_handle.DataAspectRatioMode = 'manual';
            obj.axes_handle.DataAspectRatio= [1;1;1];
            
            if par.pose_dim>3
                obj.axes_handle.CameraPosition = [-10,-22,8];
            else 
                obj.axes_handle.CameraPosition = [0,-1,0];
            end
            %obj.axes_handle.CameraTarget = [-0.913362548387552,0.616539044482742,-0.747766910811659];
            %obj.axes_handle.CameraViewAngle = 9.3176;
        end
        function obj = ResetFigureLimits(obj,limits,spacing)
            
            obj.axes_handle.XLim = limits(1,:);
            obj.axes_handle.YLim = limits(2,:);
            obj.axes_handle.ZLim = limits(3,:);
            obj.axes_handle.XTick = limits(1,1):(limits(1,2)-limits(1,1))/spacing:limits(1,2); 
            obj.axes_handle.YTick = limits(2,1):(limits(2,2)-limits(2,1))/spacing:limits(2,2); 
            obj.axes_handle.ZTick = limits(3,1):(limits(3,2)-limits(3,1))/spacing:limits(3,2);   
            
            obj.axes_handle.DataAspectRatioMode = 'manual';
            obj.axes_handle.DataAspectRatio= [1;1;1];
            
%             obj.axes_handle.CameraPosition = [-2.837916634960155,-22.911510361523757,13.494156116891816];
            obj.axes_handle.CameraPosition = [-30,-45,15];
        
        end
        function frame = SetFrame(obj,cdpr_v,cdpr_p)
            for i=1:obj.n_cables
                clearpoints(obj.lines.pulley(i));
                clearpoints(obj.lines.cable(i));
                clearpoints(obj.lines.platformPoint(i));
                clearpoints(obj.lines.platformEdge(i));
                clearpoints(obj.lines.platformPlane(i));
                clearpoints(obj.lines.PulleyAxis(i));
                R = cdpr_v.platform.rot_mat;
                k = cdpr_v.platform.rot_mat(:,3);
                r = cdpr_v.platform.position;
                d = cdpr_v.cable(i).pos_OD_glob;
                swiveling_axis = d+R*cdpr_p.cable(i).vers_k_loc*0.1;
                p(:,i) = RecordType.GetIntersectionInPlane(r,d,k);
                for j = 0:0.05:2*pi+0.05
                    vers_n_rot = R*(cdpr_v.cable(i).vers_u*cos(j) +...
                        cdpr_p.cable(i).vers_k_loc*sin(j));
                    r = d+cdpr_p.cable(i).swivel_pulley_r*...
                        (R*cdpr_v.cable(i).vers_u+vers_n_rot);
                    addpoints(obj.lines.pulley(i),[r(1)],[r(2)],[r(3)]);
                end
                addpoints(obj.lines.platformEdge(i),...
                    [d(1) cdpr_v.platform.position(1)],...
                    [d(2) cdpr_v.platform.position(2)],...
                    [d(3) cdpr_v.platform.position(3)]);
                addpoints(obj.lines.platformPoint(i),...
                    [p(1,i) cdpr_v.platform.position(1)],...
                    [p(2,i) cdpr_v.platform.position(2)],...
                    [p(3,i) cdpr_v.platform.position(3)]);
                addpoints(obj.lines.platformPlane(i),...
                    [d(1) p(1,i)],...
                    [d(2) p(2,i)],...
                    [d(3) p(3,i)]);
                r = d+cdpr_p.cable(i).swivel_pulley_r*...
                    R*(cdpr_v.cable(i).vers_u+cdpr_v.cable(i).vers_n);
                addpoints(obj.lines.cable(i),...
                    [cdpr_p.cable(i).pos_OA_glob(1) r(1,1)],...
                    [cdpr_p.cable(i).pos_OA_glob(2) r(2,1)],...
                    [cdpr_p.cable(i).pos_OA_glob(3) r(3,1)]);
                addpoints(obj.lines.PulleyAxis(i),...         % pulley axis? What is?
                    [d(1) swiveling_axis(1)],...
                    [d(2) swiveling_axis(2)],...
                    [d(3) swiveling_axis(3)]);
                frame = getframe(obj.figure_handle);
            end
        end
    end
    methods (Static)
        function p = GetIntersectionInPlane(r,d,k)
            A = zeros(3);
            v = zeros(3,1);
            A(3,:) = k';
            v(3) =  k'*d;
            [~,index] = max(abs(k));
            switch index
                case 1
                    A(1:2,:) = [k(2) -k(1) 0;
                        k(3) 0 -k(1)];
                    v(1:2,1) = [-k(1)*r(2)+k(2)*r(1);-k(1)*r(3)+k(3)*r(1)];
                case 2
                    A(1:2,:) = [k(2) -k(1) 0;
                        0 k(3) -k(2)];
                    v(1:2,1) = [-k(1)*r(2)+k(2)*r(1);-k(2)*r(3)+k(3)*r(2)];
                case 3
                    A(1:2,:) = [k(3) 0 -k(1) ;
                        0 k(3) -k(2)];
                    v(1:2,1) = [-k(1)*r(3)+k(3)*r(1);-k(2)*r(3)+k(3)*r(2)];
            end
            p = linsolve(A,v);
        end
    end
end