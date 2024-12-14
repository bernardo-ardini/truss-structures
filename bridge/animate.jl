using GLMakie
using DelimitedFiles
using Printf

e=Int.(readdlm("e.dat"))

global T=-1
global xs=[]

open("xs.dat") do f
	global T
	while !eof(f)
		s=readline(f)
		if s==@sprintf("# time %d",T+1)
			T=T+1
			open("tmp.dat","w") do file
    			write(file,"")
			end
		elseif s=="# END"
			break
		elseif s=="# end"
			push!(xs,readdlm("tmp.dat"))
		else 
			open("tmp.dat","a") do file
    			write(file,s)
    			write(file,"\n")
			end
		end
	end
end

N=1;
a=0;
b=0;

s=Scene()
x=Observable(xs[1])

for l=1:size(e,1)
	xx=@lift([$x[1,e[l,1]+1],$x[1,e[l,2]+1]])
	yy=@lift([$x[2,e[l,1]+1],$x[2,e[l,2]+1]])
	zz=@lift([$x[3,e[l,1]+1],$x[3,e[l,2]+1]])
	lines!(s,xx,yy,zz,color=:black)
end
	
record(s,"bridge.gif",0:T*N;framerate=40) do t
	x[]=xs[t%T+1]
	
	d=360
	theta=60+a*20*sin(2*pi*t/(2*T*N))
	phi=80+b*360*t/(T*N)
	
	cam=Makie.camera(s)

	w,h=size(s)
	nearplane=0.1f0
	farplane=10000f0
	aspect=Float32(w/h)
	cam.projection[]=Makie.perspectiveprojection(45f0,aspect,nearplane,farplane)
	lookat=Vec3f(0,0,20)
	eyeposition=lookat+Vec3f(d*cosd(phi)*sind(theta),d*sind(phi)*sind(theta),d*cosd(theta))
	#lookat=Vec3f(-10,0,0)
	#eyeposition=Vec3f(30,170,200)
	upvector=Vec3f(0,0,1)
	cam.view[]=Makie.lookat(eyeposition,lookat,upvector)
	
	display(t)
end