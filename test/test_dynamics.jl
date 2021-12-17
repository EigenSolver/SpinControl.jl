import Plots: scatter, plot, plot!, savefig
import Statistics: mean, std

@testset "spin dynamics" begin
    ensemble=SpinEnsemble(0.39486,3,[0,0,1],0.1,10,:spherical)
    T2=coherencetime(ensemble)
    @test abs(T2-1)<0.001
    h=50; dt=π/(h*20); 
    t=0:dt:T2

    println("generating rabi oscillation curve...")
    the_curve=analyticalrabi(t,ensemble,h)
    mc_curve=rabi(t, ensemble, h; M=400, N=500)

    error=abs.(mc_curve-the_curve)
    avg_err=mean(error); std_err=std(error);
    
    scatter(t,mc_curve,labels="Numerical", 
    fmt = :png, title="Rabi Oscillation", 
    xlabel="t",
    ylabel="f(t)")
    plot!(t, the_curve,labels="Analytical",linestyle=:dash)
    savefig("./.figs/rabi_test.png")

    println("mean error: ", avg_err)
    println("std error: ", std_err)
    @test avg_err<0.01
    @test std_err<0.01

    dt=T2/200; 
    t=0:dt:T2
    println("generating free induction decay curve...")
    the_curve=analyticalfid(t,ensemble)
    mc_curve=fid(t, ensemble; M=400, N=500)

    error=abs.(mc_curve-the_curve)
    avg_err=mean(error); std_err=std(error);
    
    scatter(t,mc_curve,labels="Numerical", 
    fmt = :png, title="Free Induction Decay", 
    xlabel="t",
    ylabel="f(t)")
    plot!(t, the_curve,labels="Analytical",linestyle=:dash)
    savefig("./.figs/fid_test.png")

    println("mean error: ", avg_err)
    println("std error: ", std_err)
    @test avg_err<0.01
    @test std_err<0.01    
end