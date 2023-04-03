%% safety critical information for modified inverted pendulum
% Or Single Integrator (added by Arthur Nascimento, June 2022)

function [h,obs_loc,circ]=obstacles_2d()
    h.nbas=1;

    [circ, obs_loc] = obstacles_1();
    x_centers = obs_loc(:,1);
    y_centers = obs_loc(:,2);
    rs = obs_loc(:,3);
    num_obs=size(circ,1);
    
    h.h=cell(num_obs,1);
    for ii=1:num_obs
        h.h{ii}= @(x) (x(1)-x_centers(ii))^2+ (x(2)-y_centers(ii))^2-rs(ii)^2;
    end
    
    h.hx=cell(num_obs,1);
    for ii=1:num_obs
        h.hx{ii}= @(x) [2*(x(1)-x_centers(ii)),2*(x(2)-y_centers(ii)),zeros(1,length(x)-2)];
    end
    
    h.hxx=cell(num_obs,1);
    for ii=1:num_obs
        h.hxx{ii}= @(x) [2 0 ;0 2];
    end
    
%     figure(1)
%     for ii=1:num_obs
%         rectangle('Position',circ(ii,:),'Curvature',[1 1],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on
%     end
end

function [circ, obs_loc] = obstacles_1()
    x_center(1) = 1.7; y_center(1) = 1.2; r(1) = 0.7;
    circ(1,:) = [x_center(1)-r(1) y_center(1)-r(1) 2*r(1) 2*r(1)];
    obs_loc=[x_center, y_center, r];
end

