% 颜色控制
Map = [1 1 1; 0 0 0];
colormap(Map);
% 设置网格大小
S = 121;
L = zeros(S);
K = eye(S);
% 把中间一个数设置为 1 作为元胞种子
M = (S+1)/2;
L(M, M) = 1;
Temp = L;
imagesc(L);
% 计算层数
Layer = (S-1)/2 + 1;
%生长率和死亡率
alpha=0.1;beta=0.05;gamma=0.05;
f=getframe(gcf);
imind=frame2im(f);
[imind,cm]=rgb2ind(imind,256);
imwrite(imind,cm,'1.gif','gif','Loopcount',inf,'DelayTime',0.1);
for t=1:Layer
    for x=M-t+1:M+t-1
          for y=M-t+1:M+t-1
              if  Temp(x,y)==1 
                  i=1; 
                  A=rand(10);
                  for m=-1:1
                      for n=-1:1
                          g=A(i);
                           %g代表growth，随机数，用来模拟繁殖概率
                  if g <= alpha && g>=0 && Temp(x+m,y+n)==0 && K(x+m,y+n)==1 
                      %向下一层成功扩散
                               Temp(x+m,y+n)=1;
                  end
                  i=i+1;
                  end
                  end
                              K(x,y)=K(x,y)-0.015;
                             if K(x,y)<=0
                                 %此处资源已被消耗完
                                 K(x,y)=0;
                               if K(x,y)==0;
                                 nd=rand(1);
                             %nd代表die of nutrition，随机数，用来模拟因资源缺乏而死亡的概率
                             if nd <= gamma && nd>=0
                                 %此处的生物死亡
                                 Temp(x,y)=0;
                             end
                               end
                             d=rand(1);
                             %d代表death，随机数，用来模拟死亡概率
                             if d <= beta && d>=0
                                 %此处的生物死亡
                                 Temp(x,y)=0;
                             end
                             end
              end
          end
    end
    L = Temp;
  imagesc(L);
    % 速度控制
   pause(0.5);
   %导出gif图
   f=getframe(gcf);
   imind=frame2im(f);
   [imind,cm]=rgb2ind(imind,256);
   imwrite(imind,cm,'1.gif','gif','WriteMode','append','DelayTime',0.1);
end
