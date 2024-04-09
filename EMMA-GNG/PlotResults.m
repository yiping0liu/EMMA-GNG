function PlotResults(S, w, C)

    N = size(w,1);
    Dim=size(S,2);
    if Dim>2
        plot3(S(:,1),S(:,2),S(:,3),'y.','MarkerSize',6);
    else
        plot(S(:,1),S(:,2),'m*');
    end
    
    hold on;
    for i=1:N-1
        for j=i:N
            if C(i,j)==1
                if Dim>2
                    plot3([w(i,1) w(j,1)],[w(i,2) w(j,2)],[w(i,3) w(j,3)],'r','LineWidth',1);
                else
                    plot([w(i,1) w(j,1)],[w(i,2) w(j,2)],'r','LineWidth',2);
                end
            end
        end
    end
    if Dim>2
          plot3(w(:,1),w(:,2),w(:,3),'o','color',[0, 0, 0.5],'MarkerFaceColor',[0.5294, 0.8078, 0.9216],'MarkerSize',6);
    else
        plot(w(:,1),w(:,2),'ko','MarkerFaceColor','y','MarkerSize',8);
    end
    hold off;  
end