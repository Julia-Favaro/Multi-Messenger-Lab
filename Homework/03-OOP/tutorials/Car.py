# Example of class describing a car
# In this case the class is defined in an external file
# you can import it from another script by using import Car from Car

class Car:
    
    #this is the constructor
    def __init__(self, m_model):
 
        self.model = m_model
        self.gear = 1    # gear is set at one
        self.speed = 0  #speed is set at zero
        
        print("%s is ready" % m_model)
        
    def speed_up(self, increase): 
        # the convention says the methods should be lower case 
        # and words can be separated by underscores
        print("the %s is speeding up!!!" % self.model)
        self.speed+=increase
        
    def check_speed(self):
        print("Speed of %s is %.f" % (self.model, self.speed))
