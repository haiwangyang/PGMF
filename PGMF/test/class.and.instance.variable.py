class Boy:
    gender = "male" # class variable
    
    def __init__(self, name):
        self.name = name # instance variable, unique to object

a = Boy("Tom")
b = Boy("Sam")

print(a.gender)
print(a.name)

print(b.gender)
print(b.name)
