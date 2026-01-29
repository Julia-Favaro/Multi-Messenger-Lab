# Multimessenger Physics Laboratory
# AA 2023-2024
# Example of class
#

class Dog():
    # this is the constructor, it is run automatically when an instance of the class is created
    # needs to be called __init__(), python recognizes it automatically as the constructor
    def __init__(
        self,  # all class method have "self" as first variable, self stands for the class itself
        name,
        breed,
        owner,
        speed=0  # default variable!
    ):
        # class attributes, valid across the whole class
        self.name = name
        self.breed = breed
        self.owner = owner
        self.speed = speed

        print(f'A dog named {self.name} is built!')
    

    def run(self, destination, speed):
        # assigning a new value to the class attribute
        self.speed = speed
        print(f'{self.name} runs towards {destination} with {self.owner} at speed {speed}!')
        
    def check_speed(self):
        print(f'{self.name} is actually running at speed {self.speed}')
    
    def visit_dog_area(self, speed=4):
        self.speed = speed
        self.run(destination='dog area', speed=self.speed)
        
        print(f'{self.name} is having fun with other dogs!')


