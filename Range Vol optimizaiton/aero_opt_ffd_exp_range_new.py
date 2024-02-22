import os
import argparse
import ast
import numpy as np
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW
from pygeo import DVGeometry, DVConstraints
from pyoptsparse import Optimization, OPT
from idwarp import USMesh
from multipoint import multiPointSparse
import pygeo
from pprint import pprint

# Use Python's built-in Argument parser to get commandline options
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, default="output")
parser.add_argument("--opt", type=str, default="IPOPT", choices=["IPOPT", "SLSQP", "SNOPT"])
parser.add_argument("--gridFile", type=str, default="wing_vol_5.cgns")
parser.add_argument("--optOptions", type=ast.literal_eval, default={}, help="additional optimizer options to be added")
args = parser.parse_args()


MP = multiPointSparse(MPI.COMM_WORLD)
MP.addProcessorSet("cruise", nMembers=1, memberSizes=MPI.COMM_WORLD.size)
comm, setComm, setFlags, groupFlags, ptID = MP.createCommunicators()

if not os.path.exists(args.output):
    if comm.rank == 0:
        os.mkdir(args.output)

aeroOptions = {
    # I/O Parameters
    "gridFile": args.gridFile,
    "outputDirectory": args.output,
    "monitorvariables": ["resrho", "cl", "cd"],
    "writeTecplotSurfaceSolution": True,
    "nSubiterTurb":8, ##CHANGED
    # Physics Parameters
    "equationType": "RANS",
    # Solver Parameters
    "smoother": "DADI",
    "MGCycle": "sg",
    "infchangecorrection": True,
    # ANK Solver Parameters
    "useANKSolver": True,
    ""
    # NK Solver Parameters
    "useNKSolver": True,
    "nkswitchtol": 1e-6,
    "ANKGlobalPreconditioner":"multigrid", ##CHANGED
    # Termination Criteria
    "L2Convergence": 1e-9,
    "L2ConvergenceCoarse": 1e-2,
    "nCycles": 10000,
    # Adjoint Parameters
    "adjointL2Convergence": 1e-10,
}

# Create solver
CFDSolver = ADFLOW(options=aeroOptions, comm=comm)
## ffd generation 
# surf = CFDSolver.getTriangulatedMeshSurface()
# # leList = [[0.01, 0, 0.001], [7.51, 0, 13.99]]
# # teList = [[4.99, 0, 0.001], [8.99, 0, 13.99]]

# # leList = [[0.0, 0.0, 0.0], [3.7538574, 0.0, 4.0572925], [11.884086, 0.0, 14.907812] ]
# # teList = [[ 8.31515, 0.02337827, 0.0], [9.2952031, -0.0096776941, 4.0572925], [14.3, 0.0, 14.907812]]
leList = [[0.1, 0.0, 0.001], [3.8, 0.0, 4.056], [11.95, 0.0, 14.90] ]
teList = [[8.31, 0.0, 0.001], [9.28, 0.0, 4.056], [14.3, 0.0, 14.90]]

# surfFormat = "point-vector"

# outFile = "wing_ffd_5_y_6.xyz"

# nSpan = [2,6]

# nChord = 8
# # relMargins = [1, 1, 1]

# # absMargins = [1, 1, 1]

# # relMargins = [0.1, 0.1, 0.1]

# # absMargins = [0.5, 0.1, 0.5]

# relMargins = [0.1, 0.05, 0.1]

# absMargins = [0.4, 0.05, 0.5]
# # relMargins = [0.01, 0.001, 0.01]

# # absMargins = [0.05, 0.001, 0.05]



# liftIndex = 2
# pygeo.geo_utils.ffd_generation.createFittedWingFFD(surf, surfFormat, outFile, leList, teList, nSpan, nChord, absMargins, relMargins, liftIndex)




CFDSolver.addLiftDistribution(150, "z")
CFDSolver.addSlices("z", np.linspace(0.1, 14, 10))


ap = AeroProblem(name="wing", alpha=0, mach=0.8, altitude=10600, areaRef=74.8, chordRef=5.7, evalFuncs=["cl", "cd"])

# Add angle of attack variable
ap.addDV("alpha", value=0.0, lower=0, upper=5.0, scale=0.1)


# Create DVGeometry object
FFDFile = "wing_ffd_5_y_6.xyz"
DVGeo = DVGeometry(FFDFile)

# Create reference axis
nRefAxPts = DVGeo.addRefAxis("wing", xFraction=0.25, alignIndex="k")
nTwist = nRefAxPts - 1


# Set up global design variables
def twist(val, geo):
    for i in range(1, nRefAxPts):
        geo.rot_z["wing"].coef[i] = val[i - 1]


DVGeo.addGlobalDV(dvName="twist", value=[0] * nTwist, func=twist, lower=-8, upper=8, scale=0.1)

# Set up local design variables
DVGeo.addLocalDV("local", lower=-0.1, upper=0.1, axis="y", scale=10)

# Add DVGeo object to CFD solver
CFDSolver.setDVGeo(DVGeo)


DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)

# Only ADflow has the getTriangulatedSurface Function
DVCon.setSurface(CFDSolver.getTriangulatedMeshSurface())

# Volume constraints
# leList = [[0.01, 0, 0.001], [7.51, 0, 13.99]]
# teList = [[4.99, 0, 0.001], [8.99, 0, 13.99]]
DVCon.addVolumeConstraint(leList, teList, nSpan=20, nChord=20, name='fuel_vol', scaled=False, addToPyOpt=False)

DVCon.addVolumeConstraint(leList, teList, nSpan=20, nChord=20, lower=1.0, scaled=True)

# Thickness constraints
DVCon.addThicknessConstraints2D(leList, teList, nSpan=10, nChord=10, lower=1.0, scaled=True)


# Le/Te constraints
DVCon.addLeTeConstraints(0, "iLow")
DVCon.addLeTeConstraints(0, "iHigh")

if comm.rank == 0:
    # Only make one processor do this
    DVCon.writeTecplot(os.path.join(args.output, "constraints.dat"))

meshOptions = {"gridFile": args.gridFile}
mesh = USMesh(options=meshOptions, comm=comm)
CFDSolver.setMesh(mesh)


def cruiseFuncs(x):
    if comm.rank == 0:
        print(x)
    # Set design vars
    DVGeo.setDesignVars(x)
    ap.setDesignVars(x)
    # Run CFD
    CFDSolver(ap)
    # Evaluate functions
    funcs = {}
    DVCon.evalFunctions(funcs)
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if comm.rank == 0:
        pprint(funcs)
    return funcs


def cruiseFuncsSens(x, funcs):
    funcsSens = {}
    DVCon.evalFunctionsSens(funcsSens)
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    CFDSolver.checkAdjointFailure(ap, funcsSens)
    return funcsSens


def objCon(funcs, printOK):

    # Wf = Wi - Wcalc
    # range = funcs[ap["cl"]]/funcs[ap["cd"]] * 2.3 * vel/9.81/SFC *ln (Wi/Wf)
    # funcs["obj"] = range
    thrust = 71200 # N
    mdot = 1.196 #kg/s 4306 kg/h
    Wi = 95000 # maximum take off mass
    rhof = 800 # kg/L
    # import pdb
    # pdb.set_trace()
    # DVCon.evalFunctions(funcs)
    # for key, value in funcs.items() :
    #     print (key, value)
    print("Velocity {}".format(ap.V))   
    #  
    # funcs["obj"] = (ap.V/9.81)*(thrust/mdot)*(funcs[ap["cl"]]/funcs[ap["cd"]])*np.log(Wi/(Wi- 2*rhof*funcs['DVCon1_volume_constraint_0']))   # Range Equation *funcs[ap["cd"]]
    
    funcs["obj"] =(-1)* (ap.V/9.81)*(thrust/mdot)*(funcs[ap["cl"]]/(funcs[ap["cd"]] + 0.04))*np.log(Wi/(Wi- 0.625*2*rhof*funcs['fuel_vol']))*(1/1000)  
    # funcs["obj"] = funcs[ap["cd"]]
    
    # Assemble the objective and any additiona constraints:
    # funcs["obj"] = funcs[ap["cd"]]/1000 ## changed
    funcs["cl_con_" + ap.name] = funcs[ap["cl"]]
    if printOK:
        print("funcs in obj:", funcs)
    return funcs


# Create optimization problem
optProb = Optimization("opt", MP.obj, comm=comm)

# Add objective
optProb.addObj("obj", scale=1e-3) #1e2

# Add variables from the AeroProblem
ap.addVariablesPyOpt(optProb)

# Add DVGeo variables
DVGeo.addVariablesPyOpt(optProb)

# Add constraints
DVCon.addConstraintsPyOpt(optProb)
optProb.addCon("cl_con_" + ap.name, lower=0.5, upper=0.5, scale=2)

# The MP object needs the 'obj' and 'sens' function for each proc set,
# the optimization problem and what the objcon function is:
MP.setProcSetObjFunc("cruise", cruiseFuncs)
MP.setProcSetSensFunc("cruise", cruiseFuncsSens)
MP.setObjCon(objCon)
MP.setOptProb(optProb)
optProb.printSparsity()



# Set up optimizer
if args.opt == "SLSQP":
    optOptions = {"IFILE": os.path.join(args.output, "SLSQP.out")}
elif args.opt == "SNOPT":
    optOptions = {
        "Major feasibility tolerance": 1e-4,
        "Major optimality tolerance": 1e-4,
        "Hessian full memory": None,
        "Function precision": 1e-8,
        "Print file": os.path.join(args.output, "SNOPT_print.out"),
        "Summary file": os.path.join(args.output, "SNOPT_summary.out"),
        "Major iterations limit": 1000,
    }
elif args.opt == "IPOPT":
    optOptions = {
        "limited_memory_max_history": 1000,
        "print_level": 5,
        "output_file": os.path.join(args.output, "IPOPT_summary.out"), # changed
        "tol": 1e-5, #1e-6,
        "acceptable_tol": 1e-5,
        "max_iter": 200, # 300
        "start_with_resto": "yes",
    }
optOptions.update(args.optOptions)
opt = OPT(args.opt, options=optOptions)

# Run Optimization
sol = opt(optProb, MP.sens, storeHistory=os.path.join(args.output, "opt.hst"))
if comm.rank == 0:
    print(sol)

