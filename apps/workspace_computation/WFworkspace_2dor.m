function out = WFworkspace_2dor(cdpr_p,cdpr_v,ut,pose,tau_lim,ws_info,out,cableIdx,varargin)
%This tension distribution algorithm was designed based on "Gouttefarde2015"

%ws= WFworkspace_2dor(Jt,S)
% %start
% ii=7;
% jj=8;
Im=[1,2,3,4,5,6,7,8];
% iii=1;
%     
% %instead of taking two random cab I switch the one with lowest cond_num
% for iii=1:7
%     ii=iii;
%     [Jc,Jd]=PermuteJ_2DOR(cdpr_v.geometric_jacobian,ii,jj);
%     tempp(iii) = cond(Jd);
% end

% [~,ii] = min(tempp);
iii=cableIdx(1);             %%% Here the choice of the force controlled cables is made before
jjj=cableIdx(2);

% [Jc,Jd]=PermuteJ_2DOR(cdpr_v.geometric_jacobian,ii,jj);
% sens.ii = ii;
% sens.jj = jj;
% sens.Jc= Jc;
% sens.Jd = Jd;

%%% FZ redefined the computation of the permutation matrix 
[Jc,Jd,P]=PermJac_2dor(cdpr_v.geometric_jacobian_l,iii,jjj);
sens.ii = iii;
sens.jj = jjj;
sens.Jc= Jc;
sens.Jd = Jd;

%%%% rows of the starting point of the algorithm must be optimized
ii=7;
jj=8;

%compute N and fp
N=[-linsolve(Jd,Jc); eye(2,2)];
cdpr_v = CalcExternalLoads(cdpr_v,cdpr_p);
fp=[linsolve(Jd,cdpr_v.platform.ext_load);0;0];
qmax=ws_info.tension_limits(2)*ones(cdpr_p.n_cables,1)-fp;
qmin=ws_info.tension_limits(1)*ones(cdpr_p.n_cables,1)-fp;

% double check for interference
% if Cab2CabInterfLumelsky(cdpr_v,cdpr_p)
    % if CableInterfPerreault(cdpr_v,cdpr_p)
    %     disp('okkk')
    % else
    %     disp('not okkk')
    % end
% end

%uncomment/comment for graphic representation
% figure()
% 
% x=linspace(-2000,15000,1000);
% plot(ws_info.tension_limits(1)*ones(1,1000),x,'k',ws_info.tension_limits(2)*ones(1,1000),x,'k')
% hold on
% plot(x,ws_info.tension_limits(1)*ones(1,1000),'k',x,ws_info.tension_limits(2)*ones(1,1000),'k')
% hold on
% for bb=1:6
%     ni=N(bb,:);
%     Limax=Lineq(x,-ni(1)/ni(2),qmax(bb)/ni(2));
%     Limin=Lineq(x,-ni(1)/ni(2),qmin(bb)/ni(2));
%     if bb==1
%         plot(x,Limax,'r',x,Limin,'r')
%         axis equal
%         hold on
%     elseif bb==2
%         plot(x,Limax,'b',x,Limin,'b')
%         hold on
%     elseif bb==3
%         plot(x,Limax,'g',x,Limin,'g')
%         hold on
%     elseif bb==4
%         plot(x,Limax,'m',x,Limin,'m')
%         hold on
%     elseif bb==5
%         plot(x,Limax,'c',x,Limin,'c')
%         hold on
%     else
%         plot(x,Limax,'y',x,Limin,'y')
%         hold on
% 
%     end
% 
% end

%determine first vertex vij intersecting two random lines Li and Lj
jj=8; %%FZ modified the initial value for the algorithm to let it converge
ni=N(ii,:);
nj=N(jj,:);

fc = linsolve([ni;nj],[qmin(ii);qmin(jj)]);
% if norm(fc)>exp(12)
%     j=5;
%     fc = linsolve([ni;nj],[qmin(i);qmin(j)]);
% else
% end

%  %not useful, underline the first lines
%             Li_min=Lineq(x,-ni(1)/ni(2),qmin(i)/ni(2));%%Lineq(x,-ni(1)/ni(2),qmin(i)/ni(2));
%             Lj_min=Lineq(x,-nj(1)/nj(2),qmin(j)/nj(2));
% %             figure(3)
%             plot(x,Li_min,'y--',x,Lj_min,'m--')
%  %uncomment/comment for graphic representation
% hold on
% plot(fc(1),fc(2),'*')
%             axis equal
%
eps1 = 10^(-8);
eps2 = 10^(-5);%-2
In=zeros(1,cdpr_p.n_cables);
for ij=1:cdpr_p.n_cables
    if ((abs(N(ij,:)*fc-qmin(ij))<eps2)||(abs(N(ij,:)*fc-qmax(ij))<eps2)||((N(ij,:)*fc-qmin(ij)>0)&&(N(ij,:)*fc-qmax(ij)<0)))%eps1
        In(ij)=ij;
    else
        In(ij)=0;
    end
end
vf=fc;
%follow line Li and compute ni_perp
aa=1;
cnt=0;
while aa

    ni=N(ii,:);
    nj=N(jj,:);
    ni_p1 = [ni(2) -ni(1)];
    ni_p2 = [-ni(2) ni(1)];
    if abs(nj*fc-qmin(jj))<eps2
        if nj*ni_p1'>=0
            ni_p=ni_p1;
        else
            ni_p=ni_p2;
        end
    elseif abs(nj*fc-qmax(jj))<eps2
        if nj*ni_p1'<=0
            ni_p=ni_p1;
        else
            ni_p=ni_p2;
        end
    else
        msgbox('miao')
    end
    alfa=zeros(1,cdpr_p.n_cables);
    for k=1:cdpr_p.n_cables
        nk=N(k,:);
        if nk*ni_p'>eps1
            if (nk*fc-qmin(k)<0 && ~(abs(nk*fc-qmin(k))<eps2))%eps1
                alfa(k)=(qmin(k)-nk*fc)/(nk*ni_p');
            elseif (((nk*fc-qmin(k)>0) && (nk*fc-qmax(k)<0)) || (abs(nk*fc-qmin(k))<eps2))%eps1
                alfa(k)=(qmax(k)-nk*fc)/(nk*ni_p');
            else
                %                 if k==ii
                %                 alfa(k)=0;
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             else
                %
                alfa(k)=10^(16);
                %             end
                %                 alfa(k)=0;%%%%%%%%%%%%%%%%%%%%%%% alfa(k)=0;
            end
        elseif ((nk*ni_p'<0)  &&  (abs(nk*ni_p')>eps1)) %%nk*ni_p'<-eps1
            if (nk*fc-qmax(k)>0 && ~(abs(nk*fc-qmax(k))<eps1))
                alfa(k)=(qmax(k)-nk*fc)/(nk*ni_p');
            elseif (((nk*fc-qmin(k)>0) && (nk*fc-qmax(k)<0))|| (abs(nk*fc-qmax(k))<eps1))
                alfa(k)=(qmin(k)-nk*fc)/(nk*ni_p');
            else
                %                 if k==ii
                %                 alfa(k)=0;
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 else
                %
                alfa(k)=10^(16);
                %                  end
                %                 alfa(k)=0;%%%%%%%%%%%%%%%%%%%%%%%%% alfa(k)=0;
            end

        else %considered lines are parallel
            if k==ii
                alfa(k)=0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else

                alfa(k)=10^(16);
            end
            %             alfa(k)=0;%%%%%%%%%%%%%%%%%%%%%%%%% alfa(k)=0;
            %             end
        end
    end
    ralfa=round(alfa,6);

    %     fi = find(ralfa==0);
    %     ff=ii;
    %     if length(fi)>1
    % %         In(fi(1))=fi(1);In(fi(2))=fi(2);
    %         nf1=N(fi(1),:);
    %         nf2=N(fi(2),:);
    %         if abs(nf1*fc-qmin(fi(1)))<eps2
    %             if nf1*ni'>0 %nf1*nj'>=0
    %                 ni=nf2;
    %                 ff=fi(2);
    %             else
    %                 ni=nf1;
    %                 ff=fi(1);
    %             end
    %         elseif abs(nf1*fc-qmax(fi(1)))<eps2
    %             if nf1*ni'<0 %nf1*nj'<=0
    %                 ni=nf2;
    %                 ff=fi(2);
    %
    %             else
    %                 ni=nf1;
    %                 ff=fi(1);
    %             end
    %         else
    %             msgbox('miao')
    %         end
    % if ff~=ii
    %         %
    %         ni_p1 = [ni(2) -ni(1)];
    %         ni_p2 = [-ni(2) ni(1)];
    %         if abs(nj*fc-qmin(jj))<eps2
    %             if nj*ni_p1'>=0
    %                 ni_p=ni_p1;
    %             else
    %                 ni_p=ni_p2;
    %             end
    %         elseif abs(nj*fc-qmax(jj))<eps2
    %             if nj*ni_p1'<=0
    %                 ni_p=ni_p1;
    %             else
    %                 ni_p=ni_p2;
    %             end
    %         else
    %             msgbox('miao')
    %         end
    %         alfa=zeros(1,3);
    %         for k=1:S.m
    %             nk=N(k,:);
    %             if nk*ni_p'>eps1
    %                 if (nk*fc-qmin(k)<0 && ~(abs(nk*fc-qmin(k))<eps1))
    %                     alfa(k)=(qmin(k)-nk*fc)/(nk*ni_p');
    %                 elseif (((nk*fc-qmin(k)>0) && (nk*fc-qmax(k)<0)) || (abs(nk*fc-qmin(k))<eps1))
    %                     alfa(k)=(qmax(k)-nk*fc)/(nk*ni_p');
    %                 else
    %                     alfa(k)=0;
    %                 end
    %             elseif nk*ni_p'<-eps1
    %                 if (nk*fc-qmax(k)>0 && ~(abs(nk*fc-qmax(k))<eps1))
    %                     alfa(k)=(qmax(k)-nk*fc)/(nk*ni_p');
    %                 elseif (((nk*fc-qmin(k)>0) && (nk*fc-qmax(k)<0))|| (abs(nk*fc-qmax(k))<eps1))
    %                     alfa(k)=(qmin(k)-nk*fc)/(nk*ni_p');
    %                 else
    %                     alfa(k)=0;
    %                 end
    %
    %             else %considered lines are parallel
    %                 if k==ii
    %                 alfa(k)=0;
    %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %             else
    %
    %                 alfa(k)=10^(16);
    %             end
    %             end
    %         end
    %         ralfa=round(alfa,2);
    %
    %         %
    %     else
    % end
    %     end
    % %     if ff==ii

    [al,~]=min(alfa(ralfa>0));
    if isempty(al)
        %                     fc=fc;
        %         i=j;
        %         j=i;
    else
        fc=fc+al*ni_p';
        ll=find(alfa==al);
    end
    if isempty(find(In==ll, ll))
        In(ll)=ll;
        vf=fc;
        jj=ii;
        ii=ll;
        aa=1;
        vertex=zeros(2,8);
    else
        vertex(:,ll)=fc;
    
        if norm(fc-vf)>eps1%%
            jj=ii;
            ii=ll;
            aa=1;

        else
            %                 if S.R==0
            %
            %                 if In(1)==Im(1) && In(2)==Im(2) && In(3)==Im(3)&& In(4)==Im(4)
            %                     %  msgbox('polygon determined')
            %                     aa=0;
            %
            %                     %                                 plot(S.r(1),S.r(2),'.','Color','b');
            %                     %                                 hold on
            %
            %                     ws=[S.r(1);S.r(2)];
            %
            %                 else
            %                     ws=[];
            %                     %                                 msgbox('polygon inexistent')
            %                     aa=0;
            %
            %                     %                                 plot(S.r(1),S.r(2),'.','Color','r');
            %                     %                                 hold on
            %
            %                 end
            %                 else
            if In(1)==Im(1) && In(2)==Im(2) && In(3)==Im(3)&& In(4)==Im(4) && In(5)==Im(5)...
                    && In(6)==Im(6) && In(7)==Im(7) && In(8)==Im(8) && ~CableInterfPerreault(cdpr_v,cdpr_p)
                %  msgbox('polygon determined')
                aa=0;
                vertex( :, all(~vertex,1) ) = [];

                out.counter = out.counter+1;
                out.pose(:,out.counter) = cdpr_v.platform.pose;
                out.position(:,out.counter) = out.pose(1:3,out.counter);
                out.ang_par(:,out.counter) = out.pose(4:end,out.counter);
                out.constr(1,out.counter) = 0;
                %%%% INTRODUCTION FZ SENSITIVITY INDEXES
                tau_c = TensDistributionMean(vertex);
                out.teiw(out.counter)=1;
                d_tau=20;       % [N]
                d_l=0.01;       % [m]
                errComb=permn([1 -1], 8)';
                % for k=1:length(errComb)
                %     err = [errComb(1:6,k).*d_l; errComb(7:8,k).*d_tau];
                %     delta_l=[err(1:6); zeros(2,1)];                               % cable length measure error
                %     delta_tau_c=err(7:8);                                         % cable tension measure error
                %     [out.sigmaTauL(out.counter),tauP_dl,~,dJ_ort]=InputRatioIndex(cdpr_v,cdpr_p,cableIdx,Jd,Jc,tau_c',delta_l);
                %     out.flag = TensionErrorInsensWS_2dor(cdpr_p,ws_info,fp,tauP_dl,N,dJ_ort,delta_tau_c,delta_l);  
                %     out.teiw(out.counter)=out.teiw(out.counter)*out.flag;          %soft-case scenario
                    out.teiw(out.counter)=norm(-Jd\Jc,inf);
%                     out.teiw(out.counter)=out.teiw(out.counter)*WrenchFeasibleErrorInsensitive(cdpr_p,ws_info,fp,tauP_dl,N,dJ_ort,delta_tau_c,delta_l);   %worst-case scenario
                    % if out.teiw(out.counter)==0
                        % break
                    % end
                % end
                [out.rotCardouSens(out.counter), out.posCardouSens(out.counter)]=SensCardou(Jd'*cdpr_v.D_mat,'inf');
                out.J_ort_cond(out.counter)=rcond(N'*N);  %only to check the condition number of the jacobian during the ws characterization
                %
            else
                %                                 msgbox('polygon inexistent')
                aa=0;
                out.flag = 0;
            end
            %                 end
        end
        %         end
        %     else
        %         aa=1;
    end

    %uncomment/comment for graphic representation
%     S.fmin = ws_info.tension_limits(1);
%     S.fmax = ws_info.tension_limits(2);
%         plot(fc(1),fc(2),'*')
%         hold on
%        ax=gca;
%         ax.XLim=[S.fmin(1)-(S.fmax(1)-S.fmin(1)) S.fmax(1)+(S.fmax(1)-S.fmin(1))];
%         ax.YLim=[S.fmin(1)-(S.fmax(1)-S.fmin(1)) S.fmax(1)+(S.fmax(1)-S.fmin(1))];
%                    axis equal
cnt=cnt+1;
    for p=1:6
        upper_lim(:,p)=linsolve([N(p,:);N(8,:)],[qmax(p);qmax(8)]);
        lower_lim(:,p)=linsolve([N(p,:);N(8,:)],[qmin(p);qmin(8)]);
    end
    delta=[max(upper_lim(1,:))-min(upper_lim(1,:)) max(lower_lim(1,:))-min(lower_lim(1,:))];
    if (cnt==1000)
        aa=0;
        out.flag = 0;
    end
end
end
% close