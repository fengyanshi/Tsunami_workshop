clear all


[x y]=meshgrid(0:20:1200-20, 0:20:1200-20);

s=10-((x-140).^2+(y-600).^2)/2.5e3;
s=max(s,0);
contourf(s)

slope=1./10.;
d=(x-100)*slope;
%d(d<=0)=0.;

save('depth.txt','d','-ascii')
save('slide.txt','s','-ascii')

