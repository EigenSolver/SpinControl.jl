@testset begin
    ensemble=SpinEnsemble(0.39486,3,[0,0,1],0.1,10,:spherical)
    @test abs(coherencetime(ensemble)-1)<0.01
    fid(ensemble)
    rabi(ensemble)
end
