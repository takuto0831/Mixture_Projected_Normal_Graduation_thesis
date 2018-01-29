a <- c(1,1,1)/sqrt(3)
cos(1/sqrt(3))
func <- function(a,b,c){
  data.frame(x=a,y=b,z=c) %>% 
    mutate(r =sqrt(x^2 + y^2 + z^2)) %>%
    mutate(x = x/r, y = y/r, z = z/r) %>% 
    mutate(theta2 = acos(z)) %>% 
    mutate(theta1 = ifelse(y/sin(theta2) > 0, acos(x/sin(theta2)), - acos(x/sin(theta2)))) %>% 
    mutate(theta1 = ifelse(theta1 < 0, theta1 + 2*pi, theta1)) %>% 
    return()
}
# real 
func(1,0,0)
func(0,0.5,2)
func(1.5,1,0)
func(1,1,1)
# pred
func(1.98,-0.15,-0.44)
func(0,0.1,0.4)
func(1.52,1.65,-0.04)
func(0.66,0.58,0.55)

# pred 2.5%
func(0.51,-1.67,-2.76)
func(-0.06,-0.04,0.33)
func(-0.33,0.64,-0.76)
func(0.40,0.41,0.42)

# 3 1 4 2
func(3.35,0.24,0.19)
func(0.05,0.23,0.54)
func(3.02,3.95,0.43)
func(0.91,0.77,0.70)


# 0.06, 0.24 ,0.06, 0.22

# 0.38, 0.30, 0.42,0.36