import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeAxisymmetricMovementProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeAxisymmetricMovementProcess(KratosMultiphysics.Process):
    """This class is used in order to impose an axisymmetric body movement in a certain region of the problem

    This class constructs the model parts containing the constrains that enforce the axisymmetric body movement
    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                        : "This process uses LinearMasterSlaveConstraint in order to impose an axisymmetric movement in the given submodelpart. The process takes the first node from the submodelpart if no node's ID is provided. The default variable is DISPLACEMENT, and in case no variable is considered for the slave the same variable will be considered",
            "computing_model_part_name"   : "computing_domain",
            "model_part_name"             : "please_specify_model_part_name",
            "new_model_part_name"         : "",
            "interval"                    : [0.0, 1e30],
            "master_variable_name"        : "DISPLACEMENT",
            "slave_variable_name"         : "",
            "axisymmetry_axis"            : [0.0,0.0,1.0],
            "max_number_of_searchs"       : 1000,
            "master_node_id"              : 0
        }
        """)

        # Detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString()):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        # The computing model part
        computing_model_part_name = settings["computing_model_part_name"].GetString()
        self.computing_model_part = Model["Structure"].GetSubModelPart(computing_model_part_name)

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # We get the corresponding model parts
        self.model_part = Model[settings["model_part_name"].GetString()]
        new_model_part_name = settings["new_model_part_name"].GetString()
        if (new_model_part_name != ""):
            if (self.model_part.HasSubModelPart(new_model_part_name)):
                self.rigid_model_part = self.model_part.GetSubModelPart(new_model_part_name)
            else:
                self.rigid_model_part = self.model_part.CreateSubModelPart(new_model_part_name)
        else:
            settings["new_model_part_name"].SetString(settings["model_part_name"].GetString())
            self.rigid_model_part = self.model_part

        # Create the process
        axisymmetric_parameters = KratosMultiphysics.Parameters("""{}""")
        axisymmetric_parameters.AddValue("model_part_name", settings["model_part_name"])
        axisymmetric_parameters.AddValue("new_model_part_name", settings["new_model_part_name"])
        axisymmetric_parameters.AddValue("master_variable_name", settings["master_variable_name"])
        axisymmetric_parameters.AddValue("slave_variable_name", settings["slave_variable_name"])
        axisymmetric_parameters.AddValue("axisymmetry_axis", settings["axisymmetry_axis"])
        axisymmetric_parameters.AddValue("max_number_of_searchs", settings["max_number_of_searchs"])
        axisymmetric_parameters.AddValue("master_node_id", settings["master_node_id"])
        self.axisymmetric_movement_process = StructuralMechanicsApplication.ImposeAxisymmetricMovementProcess(self.computing_model_part, axisymmetric_parameters)

        # Trasfering the entities
        if (new_model_part_name != ""):
            transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.rigid_model_part, self.model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.NODES)
            transfer_process.Execute()

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.axisymmetric_movement_process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        # We activate/deactivate conditions dependeding of interval
        if (self.interval.IsInInterval(current_time)):
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, self.rigid_model_part.MasterSlaveConstraints)
        else:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, False, self.rigid_model_part.MasterSlaveConstraints)
