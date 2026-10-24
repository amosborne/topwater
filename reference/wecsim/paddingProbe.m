classdef paddingProbe < bodyClass
    methods
        function obj = paddingProbe()
            h.properties = struct('centerGravity',[0 0 0],'centerBuoyancy',[0 0 0], ...
                'volume',1,'name','probe','dof',6,'dofStart',1,'dofEnd',6);
            h.simulation_parameters.direction = 0;
            h.hydro_coeffs.linear_restoring_stiffness = zeros(6);
            h.hydro_coeffs.excitation.impulse_response_fun.t = -2:2;
            h.hydro_coeffs.excitation.impulse_response_fun.f = repmat(reshape([1 2 4 8 16],1,1,5),6,1,1);
            obj@bodyClass(h);
            obj.loadHydroData(h,1);
        end

        function f = probe(obj,w)
            obj.userDefinedExcitation(w.waveAmpTime,1,0,1,1,1);
            f = obj.hydroForce.hf1.userDefinedFe(:,1);
        end

        function f = probeExtended(obj,w,dt,rampTime)
            obj.userDefinedExcitation(w.waveAmpTime,dt,w.direction,1,1,1, ...
                importdata(w.elevationFile),rampTime);
            f = obj.hydroForce.hf1.userDefinedFe(:,1);
        end

        function f = probeMultiple(obj,w,simu)
            obj.mass = 1;
            obj.hydroForcePre(w,simu,1);
            f = obj.hydroForce.hf1.userDefinedFe(:,1);
        end

        function f = probeMultipleAll(obj,w,simu)
            obj.mass = 1;
            obj.hydroForcePre(w,simu,1);
            f = obj.hydroForce.hf1.userDefinedFe;
        end

        function configure(obj,direction,t,kernel)
            h = obj.hydroData;
            h.simulation_parameters.direction = direction;
            h.hydro_coeffs.excitation.impulse_response_fun.t = t;
            h.hydro_coeffs.excitation.impulse_response_fun.f = kernel;
            obj.loadHydroData(h,1);
        end
    end

    methods (Access = protected)
        function irfInfAddedMassAndDamping(obj,varargin)
            % isolate excitation preprocessing; radiation is not under test
        end
    end
end
