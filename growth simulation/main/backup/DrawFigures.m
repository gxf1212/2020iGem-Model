function [] = DrawFigures(title,independent,dependent,label_independent,lable_dependent,microorganism,ymin)
    
    global name miu var
    
    fontweight_label="normal";
    fontweight_axis="bold";
    fontname="Arial";
    as1="B.S.";
    as2="Nostoc";
    linewidth_plot=2;
    linewidth_axis=2;
    fontsize_label=15;
    fontsize_axis=15;
    position_ylabel=[-0.1,0.5,0];
    position_xlabel=[0.5,-0.06,0];
    color_background=[1,1,1];
    if strcmp(microorganism,as1)==1
        color_curve=[1,0,0];
    end
    if strcmp(microorganism,as2)==1
        color_curve=[0,0,1];
    end

    figure('NumberTitle','off','Name',title,'color',color_background);
    
    curve_number=size(independent,1);
    if curve_number>1
        Legend(1:curve_number)="a";
        for i=1:curve_number
            if miu(i)==var
                plot(independent(i,:),dependent(i,:),'linewidth',linewidth_plot,'color',color_curve+(i-1)/curve_number*(1-color_curve));
            else
                plot(independent(i,:),dependent(i,:),'--','linewidth',linewidth_plot,'color',color_curve+(i-1)/curve_number*(1-color_curve));
            end
            
            magnitude=0;
            for magnitude=1:inf
                if miu(i)*10^magnitude>=1
                    break
                end
            end
            
            if strcmp(name,"m1")==1
                Legend(i)=strcat("m_{1}=",num2str(miu(i)*10^magnitude,'%.2f'),"×10^{",num2str(magnitude),"}");
            end
            hold on
        end
        legend(Legend,'location','northwest');
        legend('boxoff');
    end
    
    if curve_number==1
        plot(independent,dependent,'linewidth',linewidth_plot,'color',color_curve);
    end
    
    x=xlabel(label_independent,'fontname',fontname,'fontsize',fontsize_label,'fontweight',fontweight_label);
    y=ylabel(lable_dependent,'fontname',fontname,'fontsize',fontsize_label,'fontweight',fontweight_label); 
    
    set(gca,'fontsize',fontsize_axis,'linewidth',linewidth_axis,'fontweight',fontweight_axis);
    set(y,'Units','Normalized','Position',position_ylabel);
    set(x,'Units','Normalized','Position',position_xlabel);
    axis square;
    
    ylim=get(gca,'ylim');
    ylim(1)=ymin;
    set(gca,'ylim',ylim);   
    
    set(gca, 'Color', 'none');
    set(gcf,'position',[0 0 1000 1000]);
    box off;
  
end