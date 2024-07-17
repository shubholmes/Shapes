import numbers
import numpy as np
import math
from math import log as ln
import sympy as sy

def area_under_curve(function, variable, initial, final):
    """ calculates area under a curve

    Parameters
    ----------
    function : expression
        function of the curve.
    variable: a sympy.symbol object (e.g. x = sympy.symbols('x'))
        variable with respect to to be integrated
    initial: number
        start point of integration
    final: number
        end point of integration

    Returns
    -------
    no return
    """
    assert isinstance(function, sy.core.add.Add), "function must be sympy expression"
    assert isinstance(variable, sy.core.symbol.Symbol), "variable must be a sympy symbol"
    assert isinstance(initial, numbers.Number), "initial must be numeric"
    assert isinstance(final, numbers.Number), "final must be numeric"

    return float(sy.integrate(function, (variable, initial, final)).evalf())

class Rectangle:
    "Rectangle object"

    def __init__(self, length, width):
        """ instantiates the object

        Parameters
        ----------
        length: number
            length of the rectangle.
        width : number
            width of the rectangle.

        Returns
        -------
        no return
        """
        assert isinstance(length, numbers.Number), "length must be numeric"
        assert isinstance(width, numbers.Number), "width must be numeric"

        self.length = float(length)
        self.width = float(width)
        self.area = self.length * self.width
        self.perimeter = 2 * (self.length + self.width)
        self.diagonal = np.sqrt(self.length**2 + self.width**2)

    def extrude(self, height):
        assert isinstance(height, numbers.Number), "height must be numeric"
        return Cuboid(self.length, self.width, float(height))

# class Triangle_sides:
#     "Triangle object"

#     def __init__(self, sides, , base=0):
#         """ instantiates the object

#         Parameters
#         ----------
#         side_one: number
#             first side of the triangle.
#         side_two: number
#             second side of the triangle.
#         side_three: number
#             third side of the triangle.

#         Returns
#         -------
#         no return
#         """
#         assert isinstance(side_one, numbers.Number), "side_one must be numeric"
#         assert isinstance(side_two, numbers.Number), "side_two must be numeric"
#         assert isinstance(side_three, numbers.Number), "side_three must be numeric"

#         self.sides = np.array([side_one, side_two, side_three])
#         self.perimeter = sum(self.sides)
#         semi_perimeter = self.perimeter/2
#         self.area = np.sqrt(semi_perimeter*(semi_perimeter-side_one)*(semi_perimeter-side_two)*(semi_perimeter-side_three))

#         # angles are the one facing the sides
#         angle_one = math.degrees(math.acos((self.sides[1]**2 + self.sides[2]**2 - self.sides[0]**2)/(2*self.sides[1]*self.sides[2])))
#         angle_two = math.degrees(math.acos((self.sides[0]**2 + self.sides[2]**2 - self.sides[1]**2)/(2*self.sides[0]*self.sides[2])))
#         angle_three = 180 - (angle_one + angle_two)

#         self.angles = nparray([angle_one, angle_two, angle_three])
        

class Circle:
    """Circle object"""
    
    def __init__(self, diameter, points=200):
        """ instantiates the object

        Parameters
        ----------
        diameter : number
            diameter of the circle.

        Returns
        -------
        no return
        """
        assert isinstance(diameter, numbers.Number), "diameter must be numeric"
        self.diameter = float(diameter)
        self.radius = self.diameter / 2
        self.center = np.array([self.radius, self.radius])
        self.area = math.pi * (self.radius**2)
        self.perimeter = 2 * math.pi * self.radius 
        
        X = np.linspace(0, self.diameter, points)

        # r**2 = (x-h)**2 + (y-k)**2, h,k = center of circle
        # h,k = (radius, radius) for this code

        Y_pos = np.sqrt(self.radius**2 - (X-self.radius)**2) 
        Y_neg = -Y_pos
        Y_pos += self.radius
        Y_neg += self.radius

        X = np.concatenate([X[::2], np.append(X[::-2], X[0])])
        Y = np.concatenate([Y_pos[::2], np.append(Y_neg[::-2], Y_neg[0])])

        if points % 2 == 0: idx = [int(points/4)] # index of one data point to remove to make number of points equal to inputted points
        else: idx = [int(points/4), -int(points/4)] # indices of two data points to remove to make number of points equal to inputted points

        X = np.delete(X, [idx])
        Y = np.delete(Y, [idx])
        self.coordinates = np.array([X, Y])

    def arc_length(self, angle, mode="degrees"):
        """Parameters
        ----------
        angle : number
            angle of sector.
        mode : number
            mode of angle
        Returns
        -------
        arc length"""

        assert mode.casefold()=="degrees" or mode.casefold()=="radians", "mode must be either degrees or radians"
        assert isinstance(angle, numbers.Number), "angle must be numeric"

        if mode.casefold() == "radians":
            angle = math.degrees(angle)

        return ((angle/360) * self.perimeter)



    def sector(self, angle, mode="degrees"):
        """Parameters
        ----------
        angle : number
            angle of sector.
        mode : number
            mode of angle

        Returns
        -------
        area of sector
        parameter of sector"""

        assert isinstance(angle, numbers.Number), "angle must be numeric"
        assert mode.casefold()=="degrees" or mode.casefold()=="radians", "mode must be either degrees or radians"
        

        if mode.casefold() == "radians":
            angle = math.degrees(angle)

        
        area = (angle/360) * self.area
        perimeter = self.arc_length(angle) + (2 * self.radius)

        return perimeter, area

    def chord_length(self, variable_value, variable="distance", mode="degrees"):
        """Parameters
        ----------
        variable value
            distance : number
                distance of chord to the center
            angle : number
                angle of sector formed by the chord at the center of the circle.
        varaible: string
            either "distance" or "angle"
        mode : number
            mode of angle

        Returns
        -------
        chord length"""
        
        assert isinstance(variable_value, numbers.Number), "angle must be numeric"
        assert variable.casefold()=="distance" or variable.casefold()=="angle", "mode must be either degrees or radians"
        if variable == "distance":
            assert 0 < variable_value < self.radius
            length = 2 * np.sqrt(self.radius**2 - variable_value**2)
            angle = 2 * math.asin( length / (2 * self.radius) )

        else:
            assert mode.casefold()=="degrees" or mode.casefold()=="radians", "mode must be either degrees or radians"
            assert 0 < variable_value <= 360
            angle = variable_value
            if mode.casefold() == "degrees":
                angle = math.radians(angle)
            length = self.diameter * math.sin(angle/2)
        
        return math.degrees(angle), length

    def segment(self, variable_value, variable="distance", mode="degrees"):
        """Parameters
        ----------
        variable value
            distance : number
                distance of chord to the center
            angle : number
                angle of sector formed by the chord at the center of the circle.
        varaible: string
            either "distance" or "angle"
        mode : number
            mode of angle

        Returns
        -------
        area of segment
        perimeter of segment"""

        assert isinstance(variable_value, numbers.Number), "angle must be numeric"
        assert variable.casefold()=="distance" or variable.casefold()=="angle", "mode must be either degrees or radians"
        
        if variable == "distance":
            assert 0 < variable_value < self.radius

            chord_length = self.chord_length(variable_value, variable=variable, mode=mode)[1]
            angle = 2 * math.asin( chord_length / (2 * self.radius) )
           
        else:
            assert mode.casefold()=="degrees" or mode.casefold()=="radians", "mode must be either degrees or radians"
            assert 0 < variable_value <= 360

            if mode.casefold() == "degrees":
                angle = math.radians(variable_value)

            chord_length = self.chord_length(variable_value, variable=variable, mode=mode)[1]

        area = ( (self.radius**2)/2 ) * (angle - math.sin(angle))
        perimeter = chord_length + self.arc_length(angle, mode="radians") 

        return math.degrees(angle), perimeter, area

    def extrude(self, height):

        """ instantiates the object

        Parameters
        ----------
        length : number
            length of extrusion.

        Returns
        -------
        returns extruded Cylinder object
        """

        assert isinstance(length, numbers.Number), "height must be numeric"
        assert 0 < length, "length must be greater than zero"

        return Cylinder(self.diameter, height)

    def revolve(self):
        """ instantiates the object

        Returns
        -------
        returns Sphere object

        """

        assert isinstance(length, numbers.Number), "height must be numeric"
        assert 0 < length, "length must be greater than zero"
        return Sphere(self.diameter)

    def get_radius(value, parameter="area"):

        """
        Parameters
        ----------
        area : number
            area of the circle.

        Returns
        -------
        radius of the circle
        """
        assert isinstance(value, numbers.Number), "value must be numeric"
        assert parameter.casefold()=="area" or parameter.casefold()=="perimeter", "parameter must be either area or perimeter"

        if parameter.casefold() == "area":
            radius = math.sqrt( value / (math.pi) )
        else:
            radius =  value / (2 * math.pi)

        return radius



class Parabola:
    "Parabola object"

    def __init__(self, width, height, points=200, direction="down"):
        """ instantiates the object

        Parameters
        ----------
        width : number
            width of the parabola.
            Distance between the two end points of parabola
        height : number
            height of the parabola.
            Distance from base to vertex
        points : number
            number of data points for the parabola coordinates
        direction : str
            Direction parabola opens towards.
            Either "up", "down", "right" or "left"

        Returns
        -------
        no return
        """
        assert isinstance(width, numbers.Number), "width must be numeric"
        assert isinstance(height, numbers.Number), "height must be numeric"
        assert isinstance(points, int), "points must be an integer"
        assert direction.casefold() in ["up", "down", "right", "left"], "direction variable can only be 'up', 'down', 'right', 'left'."


        self.width = float(width)
        self.height = float(height)
        self.direction = direction.casefold()

        s = np.sqrt(np.square(self.width) + (16*np.square(self.height)))

        self.arc_length = ( (1/2)*s ) + ( np.square(self.width)/(8*self.height) ) * ln( ((4*self.height)+s)/self.width )


        x = sy.symbols('x')
        self.__focus_ = -self.height / (self.width- (self.width/2))**2
        expr = self.__focus_ * (x-(self.width/2))**2 + self.height
        self.area = area_under_curve(expr, x, 0, self.width)

        if self.direction == "up":
            X = np.linspace(0, self.width, points)
            self.focus = self.height / (self.width- (self.width/2))**2
            Y = self.focus * (X-(self.width/2))**2 
            self.focus_coordinates = np.array([self.width/2, self.focus])

            self.vertex = np.array([self.width/2, 0])
            self.directrix = -self.focus

        if self.direction == "down":
            X = np.linspace(0, self.width, points)
            a = -self.height / (self.width- (self.width/2))**2
            self.focus = abs(a)
            Y = -self.focus * (X-(self.width/2))**2 + self.height
            self.focus_coordinates = np.array([self.width/2, self.height-self.focus])

            self.vertex = np.array([self.width/2, self.height])
            self.directrix = self.height + self.focus

        if self.direction == "right":
            Y = np.linspace(0, self.width, (points//2)+1)
            self.focus = (self.width-(self.width/2))**2 / self.height
            X = ( (Y-(self.width/2) )**2) / self.focus 
            self.focus_coordinates = np.array([self.focus, self.width/2])

            self.vertex = np.array([0, self.width/2])
            self.directrix = -self.focus

        if self.direction == "left":
            Y = np.linspace(0, self.width, (points//2)+1)
            a = -(self.width-(self.width/2))**2 / self.height
            self.focus = abs(a)
            X = ( (Y-(self.width/2) )**2) / -self.focus + self.height
            self.focus_coordinates = np.array([self.height-self.focus, self.width/2])

            self.vertex = np.array([self.height, self.width/2])
            self.directrix = self.height + self.focus

        self.coordinates = np.array([X, Y])

    def write_script(self, name, path=r"C:\Users\Subomi\Documents\Parabola scripts\scripts"):
        path = path + '\\' + name + '.scr'
        f = open(path, "w", encoding="utf8")
        coordinates = [f"{x},{y}\n" for x, y in zip(self.coordinates[0], self.coordinates[1])]
        coordinates.insert(0, "spline\n")
        f.writelines(coordinates)
        f.close()
    
    def area_under_curve(self, initial, final):

        assert 0 <= initial <= self.width, "initial point must be within the width of the parabola"
        assert 0 <= final <= self.width, "final point must be within the width of the parabola"

        x = sy.symbols('x')
        expr = self.__focus_ * (x-(self.width/2))**2 + self.height
        integral = area_under_curve(expr, x, initial, final) # area_under_curve defined at the beginning

        return float(integral)

    def extrude(self, length):

        """ instantiates the object

        Parameters
        ----------
        length : number
            length of extrusion.

        Returns
        -------
        returns Extruded_parabola object
        """

        assert isinstance(length, numbers.Number), "height must be numeric"
        assert 0 < length, "length must be greater than zero"
        
        return Extruded_parabola(self.width, self.height, length)

class Ellipse:
    "Ellipse object"

    def __init__(self, diameters, major_axis="x", points=200):
        """ instantiates the object

        Parameters
        ----------
        major_diameter : list or arrays 
            major and minor diameters of the ellipse
        major_axis : "x", "X", "y" or "Y"
        points : number
            number of data points for the ellipse's coordinates
        "

        Returns
        -------
        no return
        """
        assert len(diameters)==2, "An ellipse has only two diameters"
        assert isinstance(diameters[0], numbers.Number), "diameter must be numeric"
        assert isinstance(diameters[1], numbers.Number), "diameter must be numeric"
        assert major_axis.isalpha(), "major_diameter must be either 'x', 'y'"
        assert major_axis.casefold() in ["x", "y"], "major_diameter must be either 'x', 'y'"
        assert isinstance(points, int), "points must be an integer"


        self.major_diameter = float(max(diameters))
        self.major_radius =  self.major_diameter/2
        self.minor_diameter = float(min(diameters))
        self.minor_radius =  self.minor_diameter/2
        self.major_axis = major_axis



        self.eccentricity = np.sqrt( 1-( self.minor_radius**2 / self.major_radius**2 ) )
        c = np.sqrt(self.major_radius**2 - self.minor_radius**2)  # c is distance between the centre on focus
        self.focus = ( self.major_diameter - ( 2*c ) ) / 2 # from major vertex


        h = (self.major_radius - self.minor_radius)**2 / (self.major_radius + self.minor_radius)**2
        self.perimeter = math.pi * ( self.major_radius + self.minor_radius ) * (1 + ( (3*h)/ ( 10 + np.sqrt(4-3*h) ) ) )
        self.area = math.pi * self.major_radius * self.minor_radius


        X = np.linspace(0, self.major_diameter, points)
        Y_pos = self.minor_radius * np.sqrt( 1 - ( (X-self.major_radius)  / self.major_radius )**2 ) # Using equation of a parabola
        Y_neg = -Y_pos

        X = np.concatenate([X[::2], np.append(X[::-2], X[0])])                         # To give a continuos circular data points
        Y = np.concatenate([Y_pos[::2], np.append(Y_neg[::-2], Y_neg[0])]) + self.minor_radius   # axes[1] added to shift y-origin to zero

        if points % 2 == 0: idx = [int(points/4)]          # index of one data point to remove to make number of points equal to inputted points
        else: idx = [int(points/4), -int(points/4)]        # indices of two data points to remove to make number of points equal to inputted points

        X = np.delete(X, [idx])
        Y = np.delete(Y, [idx])

        if self.major_axis == "x":
            self.foci_coordinates = np.array(
                    [
                        [self.focus, self.minor_radius],
                        [self.focus+(2*c), self.minor_radius]
                    ]
                )

        else:
            self.foci_coordinates = np.array(
                    [
                        [self.minor_radius, self.focus],
                        [self.minor_radius, self.focus+(2*c)]
                    ]
                )

            X, Y = Y, X

        self.coordinates = np.array([X, Y])


    def extrude(self, length):

        """ instantiates the object

        Parameters
        ----------
        length : number
            length of extrusion.

        Returns
        -------
        returns Extruded_ellipse object
        """

        assert isinstance(length, numbers.Number), "height must be numeric"
        assert 0 < length, "length must be greater than zero"
        
        return Extruded_ellipse([self.major_diameter, self.minor_diameter], length)


# Solids

class Sphere:
    "Sphere object"
        
    def __init__(self, diameter):
        """ instantiates the object

        Parameters
        ----------
        diameter : number
            diameter of the sphere.

        Returns
        -------
        no return
        """
        assert isinstance(diameter, numbers.Number), "diameter must be numeric"
        
        self.diameter = float(diameter)
        self.radius = self.diameter / 2
        self.surface_area = 4 * math.pi * self.radius**2
        self.volume = (4 * 3) * math.pi * self.radius**3

    def shell(self, shell_distance):
        """ instantiates the object

        Parameters
        ----------
        shell_distance : number
            thickness of the hollow sphere.

        Returns
        -------
        no return
        """
        assert isinstance(shell_distance, numbers.Number), "shell_distance must be numeric"
        assert self.diameter > (2*shell_distance), "Two times the shell distance must be less than the diameter of the sphere"

        self.thickness = float(shell_distance)
        self.internal_diameter = self.diameter - (2*self.thickness)
        self.internal_radius = self.internal_diameter / 2
        self.internal_surface_area = 4 * math.pi * self.internal_radius**2
        self.internal_volume = (4 * 3) * math.pi * self.internal_radius**3


class Cuboid:
    "Cube or Cuboid object"
    
    def secondary_properties(length, width, height):
        base = Rectangle(length, width)

        base_perimeter = base.perimeter
        base_area = base.area
        lateral_surface_area = base_perimeter * height
        lateral_areas = {
            "L x H": length * height,
            "W x H": width * height
        }
        total_surface_area = lateral_surface_area + (2*base_area)
        volume = base_area * height
        diagonal = np.sqrt(length**2 + width**2 + height**2)


        return base_perimeter, base_area, lateral_surface_area, lateral_areas, total_surface_area, volume, diagonal

    def __init__(self, length, width, height):
        """ instantiates the object

        Parameters
        ----------
        length :
            length of the cuboid.
        width : number
            width of the cuboid.
        height : number
            height of the cuboid.

        Returns
        -------
        no return
        """
        assert isinstance(length, numbers.Number), "length must be numeric"
        assert isinstance(width, numbers.Number), "width must be numeric"
        assert isinstance(height, numbers.Number), "height must be numeric"

        self.length = float(length)
        self.width = float(width)
        self.height = float(height)
        self.base_perimeter, self.base_area, self.lateral_surface_area, self.lateral_areas, \
        self.total_surface_area, self.volume, self.diagonal = Cuboid.secondary_properties(
                self.length, self.width, self.height
            )
        
    def shell(self, shell_distance):

        """ instantiates the object

        Parameters
        ----------
        shell_distance : number
            thickness of the hollow cuboid.

        Returns
        -------
        no return
        """

        assert isinstance(shell_distance, numbers.Number), "shell_distance must be numeric"
        assert self.length > (2*shell_distance), "Two times the shell distance must be less than the length of the cuboid"
        assert self.width > (2*shell_distance), "Two times the shell distance must be less than the width of the cuboid"
        assert self.height > (2*shell_distance), "Two times the shell distance must be less than the height of the cuboid"

        self.thickness = float(shell_distance)
        self.internal_length = self.length-(2*self.thickness)
        self.internal_width  = self.width-(2*self.thickness)
        self.internal_height = self.height-(2*self.thickness)

        self.internal_base_perimeter, self.internal_base_area, self.internal_lateral_surface_area, self.internal_lateral_areas, \
        self.internal_total_surface_area, self.internal_volume, self.internal_diagonal = Cuboid.secondary_properties(
                self.internal_length, self.internal_width, self.internal_height
            )

class Cylinder:
    "Cylinder object"
    def secondary_properties(diameter, height):
        base = Circle(diameter)
        base_perimeter = base.perimeter
        base_area = base.area
        curved_surface_area = base.perimeter * height
        total_surface_area = (2 * base.area) + curved_surface_area
        volume = base.area * height

        return base_perimeter, base_area, curved_surface_area, total_surface_area, volume

    def __init__(self, diameter, height):
        """ instantiates the object

        Parameters
        ----------
        diameter :
            diameter of the cylinder.
        height : number
            height of the cylinder.

        Returns
        -------
        no return
        """
        assert isinstance(diameter, numbers.Number), "diameter must be numeric"
        assert isinstance(height, numbers.Number), "height must be numeric"

        self.diameter = float(diameter)
        self.radius = self.diameter / 2
        self.height = float(height)

        self.base_perimeter, self.base_area, self.curved_surface_area, self.total_surface_area, \
        self.volume = Cylinder.secondary_properties(self.diameter, self.height)

    def shell(self, shell_distance):

        """ instantiates the object

        Parameters
        ----------
        shell_distance : number
            thickness of the hollow cylinder.

        Returns
        -------
        no return
        """

        assert isinstance(shell_distance, numbers.Number), "shell_distance must be numeric"
        assert self.diameter > (2*shell_distance), "Two times the shell distance must be less than the diameter of the cylinder"
        assert self.height > (2*shell_distance), "Two times the shell distance must be less than the height of the cylinder"

        self.thickness = float(shell_distance)
        self.internal_diameter = self.diameter-(2*self.thickness)
        self.internal_radius = self.internal_diameter / 2
        self.internal_height = self.height-(2*self.thickness)

        self.internal_base_perimeter, self.internal_base_area, self.internal_curved_surface_area, self.internal_total_surface_area, \
        self.internal_volume = Cylinder.secondary_properties(
                self.internal_diameter, self.internal_height
            )


class Cone:
    "Cone object"
    def secondary_properties(diameter, height, slant_height):

        base = Circle(diameter)

        base_perimeter = base.perimeter
        base_area = base.area

        curved_surface_area = math.pi * base.radius * slant_height
        total_surface_area = base.area + curved_surface_area

        volume = (1/3) * math.pi * base.radius**2 * height

        return base_perimeter, base_area, curved_surface_area, total_surface_area, volume

    def __init__(self, diameter, other_variable_value, other_variable_name="height"):
        """ instantiates the object

        Parameters
        ----------
        diameter :
            diameter of the cone.

        other_variable_value : number
            accompanying parameter value to determine other properties of the truncated cone

        other_variable_name:
            name of the accampanying parameter
            Either 'height', 'slant height' or 'angle'

        Returns
        -------
        no return
        """
        assert isinstance(diameter, numbers.Number), "diameter must be numeric"
        assert isinstance(other_variable_value, numbers.Number), "other_variable_value must be numeric"
        assert other_variable_name.casefold() in ["height", "slant_height", "angle"], "variable name must be either 'height', 'slant_height', 'angle'"

        self.diameter = float(diameter)
        self.radius = self.diameter / 2

        if other_variable_name=="height":
            self.height = float(other_variable_value)
            self.slant_height = np.sqrt(self.height**2 + self.radius**2)
            self.angle = math.degrees(math.atan(self.height/self.radius))

        elif other_variable_name=="slant_height":
            self.slant_height = float(other_variable_value)
            self.height = np.sqrt(self.slant_height**2 - self.radius**2)
            self.angle = math.degrees(math.atan(self.height/self.radius))
        else:
            assert 0<other_variable_value<90, "angle must be between 0 and 90°"
            self.angle = float(other_variable_value)
            self.height = self.radius * math.tan(math.radians(self.angle))
            self.slant_height = np.sqrt(self.height**2 + self.radius**2)


        self.base_perimeter, self.base_area, self.curved_surface_area, self.total_surface_area, self.volume = Cone.secondary_properties(
                self.diameter, self.height, self.slant_height
            )

    def frustrum(self, variable_value, variable_name="diameter"):

        assert isinstance(variable_value, numbers.Number), "variable value must be numeric"
        assert variable_name.casefold() in ["diameter", "height", "slant_height"], "variable name must be either 'diameter, height', 'slant_height'"

        if variable_name.casefold() == 'diameter':
            assert 0 < variable_value < self.diameter, "upper diameter must be less than diameter of cone."
            upper_diameter = variable_value
            upper_radius = upper_diameter/2
            height = self.height * (1 - (upper_radius/self.radius))
            return Cone_frustrum([self.diameter, upper_diameter], height, "height")

        elif variable_name.casefold() == 'height':
            assert 0 < variable_value < self.height, "height must be less than height of cone."
            upper_radius = (self.radius * (self.height - variable_value)) / self.height
            upper_diameter = 2 * upper_radius
            return Cone_frustrum([self.diameter, upper_diameter], variable_value, "height")

        else:
            assert 0 < variable_value < self.slant_height, "height must be less than slant height of cone."
            upper_diameter = 2 * (self.radius * (self.slant_height - variable_value)) / self.slant_height
            return Cone_frustrum([self.diameter, upper_diameter], variable_value, "slant_height")

    def shell(self, shell_distance):

        """ instantiates the object

        Parameters
        ----------
        shell_distance : number
            thickness of the hollow cone.

        Returns
        -------
        no return
        """
        assert isinstance(shell_distance, numbers.Number), "shell_distance must be numeric"
        assert self.diameter > (2*shell_distance), "Two times the shell distance must be less than the diameter of the cone"
        assert self.height > (2*shell_distance), "Two times the shell distance must be less than the height of the cone"
        assert self.slant_height > (2*shell_distance), "Two times the shell distance must be less than the slant_height of the cone"

        self.thickness = float(shell_distance)
        self.internal_diameter = self.diameter-(2*self.thickness)
        self.internal_radius = self.internal_diameter / 2
        self.internal_height = self.height-(2*self.thickness)

        self.internal_slant_height = np.sqrt(self.internal_height**2 + self.internal_radius**2)
        assert self.internal_slant_height > 0, "Internal slant height is less than zero"

        self.internal_base_perimeter, self.internal_base_area, self.internal_curved_surface_area, \
        self.internal_total_surface_area, self.internal_volume = Cone.secondary_properties(
                self.internal_diameter, self.internal_height, self.internal_slant_height
            )

class Cone_frustrum:
    "Cone frustrum object"
    def secondary_properties(lower_radius, upper_radius, height, slant_height, H, L):

        base = Circle(lower_radius*2)
        base_perimeter = base.perimeter
        base_area = base.area

        top = Circle(upper_radius*2)
        top_perimeter = top.perimeter
        top_area = top.area

        curved_surface_area = math.pi * ( ( lower_radius * L ) - ( upper_radius * (L-slant_height) ) )
        total_surface_area = curved_surface_area + base_area + top_area
        volume = (1/3) * math.pi * ( ( lower_radius**2 * H ) -  ( upper_radius**2 * ( H-height ) ) )

        return base_perimeter, base_area, top_perimeter, top_area, curved_surface_area, total_surface_area, volume

    def __init__(self, diameters, other_variable_value, other_variable_name="height"):

        """ instantiates the object

        Parameters
        ----------
        diameters : list or arrays 
            upper and lower diameters of the cone frustrum

        other_variable_value : number
            accompanying parameter value to determine other properties of the cone frustrum

        other_variable_name:
            name of the accampanying parameter
            Either 'height', 'slant height' or 'angle'

        Returns
        -------
        no return
        """
        assert len(diameters)==2, "Cone frustrum has only two diameters"
        assert isinstance(diameters[0], numbers.Number), "diameter must be numeric"
        assert isinstance(diameters[1], numbers.Number), "diameter must be numeric"
        assert isinstance(other_variable_value, numbers.Number), "other_variable_value must be numeric"
        assert other_variable_name.casefold() in ["height", "slant_height", "angle"], "variable name must be either 'height', 'slant_height', 'angle'"

        self.lower_diameter = float(max(diameters))
        self.lower_radius = self.lower_diameter / 2
        self.upper_diameter = float(min(diameters))
        self.upper_radius = self.upper_diameter / 2

        if other_variable_name.casefold()=="height":
            self.height = float(other_variable_value)
            self.__H = (self.lower_radius * self.height) / (self.lower_radius - self.upper_radius) #height of full cone

            self.__L = np.sqrt(self.__H**2 + self.lower_radius**2) # Slant height of full cone

            self.slant_height = self.__L - np.sqrt( (self.__H - self.height)**2 + self.upper_radius**2 )

            self.angle = math.degrees(math.asin(self.__H/self.__L))
            print(self.__H)
            print(self.__L)


        elif other_variable_name.casefold()=="slant_height":
            self.slant_height = float(other_variable_value)
            self.__L = ( self.lower_radius * self.slant_height ) / ( self.lower_radius - self.upper_radius )

            self.__H = np.sqrt(self.__L**2 - self.lower_radius**2)

            self.height = self.__H - np.sqrt( (self.__L - self.slant_height)**2 - self.upper_radius**2 ) 

            self.angle = math.degrees(math.asin(self.__H/self.__L))
            print(self.__H)
            print(self.__L)


        else:
            assert 0<other_variable_value<90, "angle must be between 0 and 90°"
            self.angle = float(other_variable_value)

            angle = math.radians(self.angle)

            self.__H = self.lower_radius * math.tan(angle)
            self.height = (self.lower_radius - self.upper_radius) * math.tan(angle)
            
            self.__L = self.lower_radius / math.cos(angle)
            self.slant_height = (self.lower_radius - self.upper_radius) / math.cos(angle)
            print(self.__H)
            print(self.__L)


        self.base_perimeter, self.base_area, self.top_perimeter, self.top_area, self.curved_surface_area, \
        self.total_surface_area, self.volume = Cone_frustrum.secondary_properties(
                self.lower_radius, self.upper_radius, self.height, self.slant_height, self.__H, self.__L
            )

    def shell(self, shell_distance):

        """ instantiates the object

        Parameters
        ----------
        shell_distance : number
            thickness of the hollow cone frustrum.

        Returns
        -------
        no return
        """

        assert isinstance(shell_distance, numbers.Number), "shell_distance must be numeric"
        assert self.lower_diameter > (2*shell_distance), "Two times the shell distance must be less than the lower diameter of the cone frustrum"
        assert self.upper_diameter > (2*shell_distance), "Two times the shell distance must be less than the upper diameter of the cone frustrum"
        assert self.height > (2*shell_distance), "Two times the shell distance must be less than the height of the cone frustrum"
        assert self.slant_height > (2*shell_distance), "Two times the shell distance must be less than the slant height of the cone frustrum"

        self.thickness = float(shell_distance)
        self.internal_lower_diameter = self.lower_diameter-(2*self.thickness)
        self.internal_lower_radius = self.internal_lower_diameter / 2
        self.internal_upper_diameter = self.upper_diameter-(2*self.thickness)
        self.internal_upper_radius = self.internal_upper_diameter / 2
        self.internal_height = self.height-(2*self.thickness)
        self.internal_slant_height = self.slant_height-(2*self.thickness)
        self.__internal_H = self.__H - (2*self.thickness)
        self.__internal_L = self.__L - (2*self.thickness)

        self.internal_base_perimeter, self.internal_base_area, self.internal_top_perimeter, self.internal_top_area, self.internal_curved_surface_area, \
        self.internal_total_surface_area, self.internal_volume = Cone_frustrum.secondary_properties(
                self.internal_lower_radius, self.internal_upper_radius, self.internal_height, 
                self.internal_slant_height, self.__internal_H, self.__internal_L
            )

class Extruded_ellipse():

    "Extruded Parabola object"
    def secondary_properties(diameters, length):

        face = Ellipse(diameters)
        curved_surface = Rectangle(length, face.perimeter)

        face_perimeter = face.perimeter
        face_area = face.area

        curved_suface_area = curved_surface.area

        total_surface_area = (2*face_area) + curved_suface_area

        volume = face_area * length

        return face_perimeter, face_area, curved_suface_area,total_surface_area, volume


    def __init__(self, diameters, length):

        """ instantiates the object

        Parameters
        ----------
        diameters : list or arrays 
            major and minor diameters of the ellipse
        major_axis : "x", "X", "y" or "Y"
        points : number
            number of data points for the ellipse's coordinates
        length : number
            length of extruded parabola

        Returns
        -------
        no return
        """

        assert len(diameters)==2, "An ellipse has only two diameters"
        assert isinstance(diameters[0], numbers.Number), "diameter must be numeric"
        assert isinstance(diameters[1], numbers.Number), "diameter must be numeric"
        assert isinstance(length, numbers.Number), "length must be numeric"
        
        self.major_diameter = float(max(diameters))
        self.major_radius =  self.major_diameter/2
        self.minor_diameter = float(min(diameters))
        self.minor_radius =  self.minor_diameter/2
        self.length = float(length)

        self.face_perimeter, self.face_area, self.curved_suface_area, self.total_surface_area, self.volume = Extruded_ellipse.secondary_properties(
                [self.major_diameter, self.minor_diameter], self.length
            )

    def shell(self, shell_distance):

        """ instantiates the object

        Parameters
        ----------
        shell_distance : number
            thickness of the hollow extruded ellipse.

        Returns
        -------
        no return
        """

        assert isinstance(shell_distance, numbers.Number), "shell_distance must be numeric"
        assert self.major_diameter > (2*shell_distance), "Two times the shell distance must be less than the major_diameter of the extruded ellipse"
        assert self.minor_diameter > (2*shell_distance), "Two times the shell distance must be less than the minor_diameter of the extruded ellipse"
        assert self.length > (2*shell_distance), "Two times the shell distance must be less than the length of the extruded ellipse"

        self.thickness = float(shell_distance)
        self.internal_major_diameter = self.major_diameter-(2*self.thickness)
        self.internal_minor_diameter = self.minor_diameter-(2*self.thickness)
        self.internal_length = self.length-(2*self.thickness)

        self.internal_face_perimeter, self.internal_face_area, self.internal_curved_suface_area, self.internal_total_surface_area, \
        self.internal_volume = Extruded_ellipse.secondary_properties(
                [self.internal_major_diameter, self.internal_minor_diameter], self.length
            )

class Extruded_parabola():
    "Extruded Parabola object"

    def secondary_properties(width, height, length):
        face = Parabola(width, height)
        curved_surface = Rectangle(length, face.arc_length)
        base = Rectangle(length, width)

        base_perimeter = base.perimeter
        curved_suface_perimeter = curved_surface.perimeter
        vertical_surface_perimeter = face.arc_length + face.width

        base_area = base.area
        curved_suface_area = curved_surface.area
        vertical_surface_area = face.area
        total_surface_area = base_area + curved_suface_area + (2*vertical_surface_area)

        volume = vertical_surface_area * length

        return base_perimeter, curved_suface_perimeter, vertical_surface_perimeter, base_area, \
        curved_suface_area, vertical_surface_area, total_surface_area, volume


    def __init__(self, width, height, length):

        """ instantiates the object

        Parameters
        ----------
        width : number
            width of the parabola.
            Distance between the two end points of parabola
        height : number
            height of the parabola.
            Distance from base to vertex
        length : number
            length of extruded parabola

        Returns
        -------
        no return
        """
        assert isinstance(width, numbers.Number), "width must be numeric"
        assert isinstance(height, numbers.Number), "height must be numeric"
        assert isinstance(length, numbers.Number), "length must be numeric"
        
        self.width = float(width)
        self.height = float(height)
        self.length = float(length)

        self.base_perimeter, self.curved_suface_perimeter, self.vertical_surface_perimeter, self.base_area, \
        self.curved_suface_area, self.vertical_surface_area, self.total_surface_area, self.volume = Extruded_parabola.secondary_properties(
                self.width, self.height, self.length
            )

    def shell(self, shell_distance):

        """ instantiates the object

        Parameters
        ----------
        shell_distance : number
            thickness of the hollow extruded parabola.

        Returns
        -------
        no return
        """

        assert isinstance(shell_distance, numbers.Number), "shell_distance must be numeric"
        assert self.width > (2*shell_distance), "Two times the shell distance must be less than the width of the extruded parabola"
        assert self.height > (2*shell_distance), "Two times the shell distance must be less than the height of the extruded parabola"
        assert self.length > (2*shell_distance), "Two times the shell distance must be less than the length of the extruded parabola"

        self.thickness = float(shell_distance)
        self.internal_width  = self.width-(2*self.thickness)
        self.internal_height = self.height-(2*self.thickness)
        self.internal_length = self.length-(2*self.thickness)

        self.internal_base_perimeter, self.internal_curved_suface_perimeter, self.internal_vertical_surface_perimeter, \
        self.internal_base_area, self.internal_curved_suface_area,  self.internal_vertical_surface_area, self.internal_total_surface_area, \
        self.internal_volume = Extruded_parabola.secondary_properties(
                self.internal_width, self.internal_height, self.internal_length
            )
