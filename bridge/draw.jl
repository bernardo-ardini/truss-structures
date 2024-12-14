using CairoMakie
using DelimitedFiles
using Printf

d=450
theta=60
phi=70

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

for i=1:T+1

s=Scene(size=(600,600))
cam=Makie.camera(s)

w,h=size(s)
nearplane=0.1f0
farplane=10000f0
aspect=Float32(w/h)
cam.projection[]=Makie.perspectiveprojection(45f0,aspect,nearplane,farplane)
lookat=Vec3f(0,0,20)
eyeposition=lookat+Vec3f(d*cosd(phi)*sind(theta),d*sind(phi)*sind(theta),d*cosd(theta))
upvector=Vec3f(0,0,1)
cam.view[]=Makie.lookat(eyeposition,lookat,upvector)

x=xs[i]

for l=1:size(e,1)
	x1=x[1,e[l,1]+1]
	x2=x[1,e[l,2]+1]
	y1=x[2,e[l,1]+1]
	y2=x[2,e[l,2]+1]
	z1=x[3,e[l,1]+1]
	z2=x[3,e[l,2]+1]
	lines!(s,[x1,x2],[y1,y2],[z1,z2],color=:black)
end

save(@sprintf("../presentation/immagini/bridge%d.pdf",i),s)

end