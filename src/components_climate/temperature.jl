# -----------------------------------------------------------
# Global Surface Temperature Change
# -----------------------------------------------------------

# ============================================================================
# 【temperature 组件】FaIR v2.0 三热箱能量平衡模块：把总辐射强迫 F 转为
# 全球平均地表温度距平 T。三个热箱分别代表不同响应时间尺度（深海平衡/
# 上层海洋调整等）：Tj[t] = Tj[t-1]·decay + F·q·(1−decay)，
# T 取相邻两期热箱和的平均。T 是下游 MESMER 降尺度与损害模块的总驱动。
# ============================================================================
@defcomp temperature begin

    #d   = Parameter(index=[3])      # Thermal response timescales: [1] Thermal equilibration of deep ocean & [2] Thermal admustment of upper ocean (years).
    decay_factor = Parameter(index=[3]) # Thermal response decay factor, calculated as exp(-1/d) where d represents thermal response timescale of jth thermal box.
    q       = Parameter(index=[3])      # Raditive forcing coefficient: [1] Thermal equilibration of deep ocean & [2] Thermal admustment of upper ocean (K W⁻¹m²).
    F       = Parameter(index=[time])   # Total radiative forcing (Wm⁻²).
    Tj_0    = Parameter(index=[3])
    T_0     = Parameter()

    T   = Variable(index=[time])    # Global mean surface temperature anomaly (K).
    Tj  = Variable(index=[time,3])  # Temperature change for three thermal pools (K).

    function run_timestep(p, v, d, t)

        if is_first(t)

        	# Set initial condition for three thermal boxes.
        	v.Tj[t,:] = p.Tj_0

        	# Set initial temperature.
        	v.T[t] = p.T_0

        else

        	#Calculate temperature change for the three different thermal response times.
            for j=1:3
                v.Tj[t,j] = v.Tj[t-1,j] * p.decay_factor[j] + p.F[t] * p.q[j] * (1.0 - p.decay_factor[j])
            end

            #Calculate global mean surface temperature anomaly.
            v.T[t] = sum([v.Tj[t-1,:] v.Tj[t,:]]) / 2
        end
    end
end
