#Example of class describing a bicycle
# In this case the class is defined in an external file

class Bicycle:
    
    #this is the constructor
    def __init__(self,m_gear,m_speed=0): #speed is given by default
        self.gear = m_gear
        self.speed = m_speed
        
    def speed_up(self,increase): 
        # the convention says the methods should be lower case 
        # and words can be separated by underscores
        print("Speeding up!!!")
        self.speed+=increase
        
    def check_speed(self):
        print("Speed is %.f" % self.speed)
