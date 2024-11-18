
! 物理定数

module physical_constants
    
    implicit none 
  
    double precision, parameter :: PI                        = 4.0d0 * atan(1.0d0)        ! 円周率
    double precision, parameter :: LIGHT_SPEED               = 2.99792458d10              ! 光速
    double precision, parameter :: GRAVITATIONAL_CONSTANT    = 6.6743d-8                  ! 万有引力定数
    double precision, parameter :: PLANCK_CONSTANT           = 6.6260755d-27              ! プランク定数
    double precision, parameter :: DIRAC_CONSTANT            = PLANCK_CONSTANT/(2 * PI)   ! ディラック定数
    double precision, parameter :: BOLTZMANN_CONSTANT        = 1.380658d-16               ! ボルツマン定数
    double precision, parameter :: ATOMIC_MASS_UNIT          = 1.6605d-24                 ! 原子質量単位
    double precision, parameter :: PROTON_MASS               = 1.6726d-24                 ! 陽子質量
    double precision, parameter :: NEUTRON_MASS              = 1.6749d-14                 ! 中性子質量
    double precision, parameter :: ELECTRON_MASS             = 9.1093897d-28              ! 電子質量
    double precision, parameter :: HYDROGEN_MASS             = 1.67353d-24                ! 水素質量
    double precision, parameter :: CHARGE_UNIT               = 4.803204673d-10            ! 電荷単位
    double precision, parameter :: ELECTRON_VOLT             = 1.60218d-12                ! 1eV
    double precision, parameter :: BOHR_RADIUS               = 5.2918d-9                  ! ボーア半径
    double precision, parameter :: CLASSICAL_ELECTRON_RADIUS = 2.8179d-13                 ! 古典電子半径
    double precision, parameter :: AVOGADRO_CONSTANT         = 6.0221d23                  ! アボガドロ定数
    double precision, parameter :: GAS_CONSTANT_MOL          = 8.3145d7                   ! 1molの気体定数
    double precision, parameter :: STEFAN_BOLTZMANN_CONSTANT = 5.6705d-5                  ! ステファン-ボルツマン定数
    double precision, parameter :: RADIATION_CONSTANT        = 7.5646d-15                 ! 輻射定数
    double precision, parameter :: ASTRONOMICAL_UNIT         = 1.496d+13                  ! 天文単位(cm)
    double precision, parameter :: LIGHT_YEAR                = 9.4605d17                  ! 光年(cm)
    double precision, parameter :: PARSEC                    = 3.0857d18                  ! パーセク(cm)
    double precision, parameter :: SOLAR_YEAR                = 3.1557d7                   ! 太陽年(s)
    double precision, parameter :: SOLAR_MASS                = 1.989d33                   ! 太陽質量(g)
    double precision, parameter :: SOLAR_RADIUS              = 6.960d10                   ! 太陽半径(cm)
    double precision, parameter :: POTASSIUM_MASS            = 39d0*PROTON_MASS             ! カリウム質量
    double precision, parameter :: K_IONIZATION_POTENTIAL    = 4.34d0 * ELECTRON_VOLT       ! カリウムのイオン化エネルギー

end module Physical_constants