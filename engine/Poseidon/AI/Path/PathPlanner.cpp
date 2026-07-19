#include <Poseidon/Core/Application.hpp>

#include <Poseidon/AI/Path/PathPlanner.hpp>
#include <Poseidon/AI/VehicleAI.hpp>
#include <Poseidon/World/World.hpp>
#include <Poseidon/World/Terrain/Landscape.hpp>
#include <Poseidon/AI/AI.hpp>

#include <Poseidon/AI/Path/AIDefs.hpp>

#include <Poseidon/World/Terrain/Roads.hpp>
#include <float.h>
#include <stdlib.h>
#include <cmath>
#include <Poseidon/Foundation/Common/FltOpts.hpp>
#include <Poseidon/Foundation/Containers/Array.hpp>
#include <Poseidon/Foundation/Containers/HashMap.hpp>
#include <Poseidon/Foundation/Enums/EnumNames.hpp>
#include <Poseidon/Foundation/Framework/DebugLog.hpp>
#include <Poseidon/Foundation/Math/MathDefs.hpp>
#include <Poseidon/Foundation/Math/MathOpt.hpp>
#include <Poseidon/Foundation/Memory/MemAlloc.hpp>
#include <Poseidon/Foundation/Types/LLinks.hpp>

#include <Poseidon/AI/Path/AStar.hpp>

namespace Poseidon
{
using namespace Foundation;
using Foundation::EnumName;

// A* algoritm for searching the best path
#define ITER_PER_CYCLE 500 // per-think budget; 50 was tuned for 2001-era CPUs
#define MAX_ITER 10000
#define Directions 20
#define direction_delta directions20

#define CHECK_PERFORMANCE 0
#if CHECK_PERFORMANCE
#undef LOG_STRAT
#define LOG_STRAT 1
#endif

union ASSField
{
    int key;
    struct
    {
        short int x;
        short int z;
    } coord;
    ASSField(int x, int z)
    {
        coord.x = x;
        coord.z = z;
    }
    bool operator==(const ASSField& with) const { return with.key == key; }
    int GetKey() const { return key; }
};

typedef AStarNode<ASSField> ASSNode;
template <>
DEFINE_FAST_ALLOCATOR(ASSNode)

class ASSCostFunction
{
  protected:
    CostFunction _costFunction;
    void* _param;
    OLink<EntityAI> _vehicle;

  public:
    ASSCostFunction(CostFunction costFunction, void* param, EntityAI* vehicle)
    {
        _costFunction = costFunction;
        _param = param;
        _vehicle = vehicle;
    }
    float operator()(const ASSField& field1, const ASSField& field2) const;

  protected:
    float GetFieldCost(int x, int z) const { return _costFunction(x, z, _param); }
};

float ASSCostFunction::operator()(const ASSField& field1, const ASSField& field2) const
{
    const float H_SQRT5 = 2.2360679775;
    const float H_SQRT5_4 = 0.25 * H_SQRT5;
    const float H_SQRT2_2 = 0.5 * H_SQRT2;
    const float ROAD_BONUS = 100.0F;

    int xf = field1.coord.x;
    int zf = field1.coord.z;
    int xt = field2.coord.x;
    int zt = field2.coord.z;

    int dx = xt - xf;
    int dz = zt - zf;

    float result = GetFieldCost(xf, zf) + GetFieldCost(xt, zt);
    switch (dx)
    {
        case -2:
            switch (dz)
            {
                case -1:
                    result += GetFieldCost(xf - 1, zf - 1);
                    result += GetFieldCost(xf - 1, zf);
                    result *= H_SQRT5_4;
                    break;
                case 0:
                    result *= 0.5;
                    result += GetFieldCost(xf - 1, zf);
                    break;
                case 1:
                    result += GetFieldCost(xf - 1, zf);
                    result += GetFieldCost(xf - 1, zf + 1);
                    result *= H_SQRT5_4;
                    break;
                case -2:
                case 2:
                    Fail("Unused");
                    break;
            }
            break;
        case -1:
            switch (dz)
            {
                case -2:
                    result += GetFieldCost(xf - 1, zf - 1);
                    result += GetFieldCost(xf, zf - 1);
                    result *= H_SQRT5_4;
                    break;
                case -1:
                    result *= H_SQRT2_2;
                    break;
                case 0:
                    result *= 0.5;
                    break;
                case 1:
                    result *= H_SQRT2_2;
                    break;
                case 2:
                    result += GetFieldCost(xf - 1, zf + 1);
                    result += GetFieldCost(xf, zf + 1);
                    result *= H_SQRT5_4;
                    break;
            }
            break;
        case 0:
            switch (dz)
            {
                case -2:
                    result *= 0.5;
                    result += GetFieldCost(xf, zf - 1);
                    break;
                case -1:
                    result *= 0.5;
                    break;
                case 1:
                    result *= 0.5;
                    break;
                case 2:
                    result *= 0.5;
                    result += GetFieldCost(xf, zf + 1);
                    break;
                case 0:
                    Fail("Unused");
                    break;
            }
            break;
        case 1:
            switch (dz)
            {
                case -2:
                    result += GetFieldCost(xf, zf - 1);
                    result += GetFieldCost(xf + 1, zf - 1);
                    result *= H_SQRT5_4;
                    break;
                case -1:
                    result *= H_SQRT2_2;
                    break;
                case 0:
                    result *= 0.5;
                    break;
                case 1:
                    result *= H_SQRT2_2;
                    break;
                case 2:
                    result += GetFieldCost(xf, zf + 1);
                    result += GetFieldCost(xf + 1, zf + 1);
                    result *= H_SQRT5_4;
                    break;
            }
            break;
        case 2:
            switch (dz)
            {
                case -1:
                    result += GetFieldCost(xf + 1, zf - 1);
                    result += GetFieldCost(xf + 1, zf);
                    result *= H_SQRT5_4;
                    break;
                case 0:
                    result *= 0.5;
                    result += GetFieldCost(xf + 1, zf);
                    break;
                case 1:
                    result += GetFieldCost(xf + 1, zf);
                    result += GetFieldCost(xf + 1, zf + 1);
                    result *= H_SQRT5_4;
                    break;
                case -2:
                case 2:
                    Fail("Unused");
                    break;
            }
            break;
    }

    EntityAI* veh = _vehicle;
    if (result < GET_UNACCESSIBLE && !veh->GetType()->IsKindOf(GLOB_WORLD->Preloaded(VTypeAir)) &&
        !veh->GetType()->IsKindOf(GLOB_WORLD->Preloaded(VTypeShip)) && !veh->IsCautious())
    {
        GeographyInfo gf = GLOB_LAND->GetGeography(xf, zf);
        GeographyInfo gt = GLOB_LAND->GetGeography(xt, zt);
        if (gf.u.road)
        {
            if (!gt.u.road)
            {
                // was on road, will be out of road
                result += ROAD_BONUS;
            }
        }
        else // !gf.road
        {
            if (gt.u.road)
            {
                // was out of road, will be on road
                result -= ROAD_BONUS;
            }
        }
        // the bonus must not make the edge negative - A* assumes non-negative costs
        saturateMax(result, 0.0F);
    }
    return result;
}

class ASSHeuristicFunction
{
  protected:
    float _coef;

  public:
    ASSHeuristicFunction(float coef) { _coef = coef; }
    float operator()(const ASSField& field1, const ASSField& field2) const
    {
        int dx = abs((int)field1.coord.x - field2.coord.x);
        int dz = abs((int)field1.coord.z - field2.coord.z);
        /*
                float minD = (dx + dz - fabs(dx - dz)) * 0.5;
                float maxD = (dx + dz + fabs(dx - dz)) * 0.5;
                return _coef * ((maxD - minD) + H_SQRT2 * minD);
        */
        // optimization
        float dif = fabs(dx - dz);
        return _coef * (dif + (0.5f * H_SQRT2) * (dx + dz - dif));
    }
};

class ASSIterator
{
  protected:
    int _x, _z;
    int _index;

  public:
    ASSIterator(const ASSField& field, void* context = nullptr)
    {
        _index = 0;
        _x = field.coord.x;
        _z = field.coord.z;
    }
    operator bool() const { return _index < Directions; }
    void operator++() { _index++; }
    operator ASSField() { return ASSField(_x + direction_delta[_index][0], _z + direction_delta[_index][1]); }
};

struct ASSOpenListTraits
{
    static bool IsLess(const ASSNode* a, const ASSNode* b) { return a->_f < b->_f; }
    static bool IsLessOrEqual(const ASSNode* a, const ASSNode* b) { return a->_f <= b->_f; }
};
class ASSOpenList : public HeapArray<ASSNode*, MemAllocD, ASSOpenListTraits>
{
    typedef HeapArray<ASSNode*, MemAllocD, ASSOpenListTraits> base;

  public:
    void UpdateUp(ASSNode* node) { base::HeapUpdateUp(node); }
    void Add(ASSNode* node) { base::HeapInsert(node); }
    bool RemoveFirst(ASSNode*& node) { return base::HeapRemoveFirst(node); }
    const ASSNode* GetFirst() const { return this->Size() > 0 ? (*this)[0] : nullptr; }
};

class ASSNodeRef : public SRef<ASSNode>
{
    typedef SRef<ASSNode> base;

  public:
    ASSNodeRef() = default;
    ASSNodeRef(ASSNode* node) : base(node) {}
    int GetKey() const { return (*this)->_field.key; }
};

struct ASSClosedListTraits
{
    typedef int KeyType;
    static unsigned int CalculateHashValue(KeyType key) { return (unsigned int)key; }

    // return negative when k1<k2, positive when k1>k2, zero when equal
    static int CmpKey(KeyType k1, KeyType k2) { return k1 - k2; }
    static KeyType GetKey(const ASSNodeRef& item) { return item.GetKey(); }
};
class ASSClosedList : public MapStringToClass<ASSNodeRef, AutoArray<ASSNodeRef>, ASSClosedListTraits>
{
    typedef MapStringToClass<ASSNodeRef, AutoArray<ASSNodeRef>, ASSClosedListTraits> base;

  public:
    int Add(ASSNode* node) { return base::Add(node); }
};

typedef AStar<ASSField, ASSCostFunction, ASSHeuristicFunction, ASSIterator, ASSClosedList, ASSOpenList> AStarStrategic;

class AIPathPlanner : public IAIPathPlanner
{
  protected:
    // strategic plan
    AutoArray<FieldPassing> _plan;
    AutoArray<Vector3> _planPoints;

    // ultimate destination
    Vector3 _destination;

    // implementation of A*
    SRef<AStarStrategic> _algorithm;
    int _iterTotal;

    float _heuristic;
    bool _searching;

    CostFunction _costFunction;
    void* _param;
    OLink<VehicleWithAI> _vehicle;

#if CHECK_PERFORMANCE
    __int64 _perfStart;
    int _timeTotal;
#endif

  public:
    AIPathPlanner(CostFunction func, void* param);
    void Init() override;
    LSError Serialize(ParamArchive& ar) override;

    bool IsSearching() const override { return _searching; }

    int GetPlanSize() const override { return _plan.Size(); }
    float GetTotalCost() const override;

    int FindBestIndex(Vector3Par pos) const override;
    bool GetPlanPosition(int index, Vector3& pos) const override;
    FieldPassing::Mode GetPlanMode(int index) const override;
    GeographyInfo GetGeography(int index) const override;
    bool IsOnPath(int x, int z, int from, int to) const override;

    bool StartSearching(AI::ThinkImportance prec, VehicleWithAI* veh, Vector3Par ptStart, Vector3Par ptEnd) override;
    void StopSearching() override;
    bool ProcessSearching() override;

  protected:
    float CalculateHeuristic(int xs, int zs, int xe, int ze);
    float GetFieldCost(int x, int z) { return _costFunction(x, z, _param); }
    bool FindNearestSafe(int& x, int& z, float threshold);
    void CalculatePlanPositions();

    float GetExposure(int x, int z);
    bool IsSafe(int x, int z, float threshold)
    {
        return InRange(x, z) && GetFieldCost(x, z) < GET_UNACCESSIBLE && GetExposure(x, z) < threshold;
    }

    float Distance(int index, Vector3Par pos) const;
};

IAIPathPlanner* CreateAIPathPlanner(CostFunction func, void* param)
{
    return new AIPathPlanner(func, param);
}

AIPathPlanner::AIPathPlanner(CostFunction func, void* param)
{

    _costFunction = func;
    _param = param;

    _iterTotal = 0;

    Init();
}

void AIPathPlanner::Init()
{
    _plan.Clear();
    _planPoints.Clear();
    _destination = VZero;
    _searching = false;
}

template <>
const EnumName* Foundation::GetEnumNames(FieldPassing::Mode dummy)
{
    static const EnumName FPModeNames[] = {EnumName(FieldPassing::Move, "MOVE"),
                                           EnumName(FieldPassing::MoveOnRoad, "ROAD"), EnumName()};
    return FPModeNames;
}

LSError FieldPassing::Serialize(ParamArchive& ar)
{
    PARAM_CHECK(ar.Serialize("x", _x, 1))
    PARAM_CHECK(ar.Serialize("z", _z, 1))
    PARAM_CHECK(ar.SerializeEnum("mode", _mode, 1, Move))
    PARAM_CHECK(ar.Serialize("cost", _cost, 1))
    return LSOK;
}

LSError AIPathPlanner::Serialize(ParamArchive& ar)
{
    PARAM_CHECK(ar.Serialize("Plan", _plan, 1))
    PARAM_CHECK(ar.SerializeArray("Points", _planPoints, 1))
    PARAM_CHECK(ar.Serialize("destination", _destination, 1))
    if (ar.IsLoading() && ar.GetPass() == ParamArchive::PassFirst)
    {
        // do not save/load partial searching
        _searching = false;
        // _costFunction, _param already initialized
    }
    return LSOK;
}

bool AIPathPlanner::StartSearching(AI::ThinkImportance prec, VehicleWithAI* veh, Vector3Par ptStart, Vector3Par ptEnd)
{
    if (!veh)
    {
        _searching = false;
        return false;
    }
    _vehicle = veh;
    _destination = ptEnd;

    MapCoord xs = toIntFloor(ptStart.X() * InvLandGrid);
    MapCoord zs = toIntFloor(ptStart.Z() * InvLandGrid);
    MapCoord xe = toIntFloor(ptEnd.X() * InvLandGrid);
    MapCoord ze = toIntFloor(ptEnd.Z() * InvLandGrid);

    saturate(xs, 0, LandRange - 1);
    saturate(zs, 0, LandRange - 1);
    saturate(xe, 0, LandRange - 1);
    saturate(ze, 0, LandRange - 1);

#if LOG_STRAT
    LOG_DEBUG(AI, "Start searching: {}, from {:.0f}, {:.0f} to {:.0f}, {:.0f}", (const char*)_vehicle->GetDebugName(),
              xs * LandGrid + 0.5 * LandGrid, zs * LandGrid + 0.5 * LandGrid, xe * LandGrid + 0.5 * LandGrid,
              ze * LandGrid + 0.5 * LandGrid);
#endif

    AI_ERROR(!_searching);
    if (prec == AI::LevelCommands ||
        (ptStart.Distance2(ptEnd) < Square(DIST_MIN_OPER) && !veh->GetType()->IsKindOf(GWorld->Preloaded(VTypeAir)) &&
         !veh->GetType()->IsKindOf(GWorld->Preloaded(VTypeShip))))
    {
        // no path will be searched
        _plan.Resize(2);
        _plan[0]._x = xs;
        _plan[0]._z = zs;
        _plan[0]._mode = FieldPassing::Move;
        _plan[1]._x = xe;
        _plan[1]._z = ze;
        _plan[1]._mode = FieldPassing::Move;
        CalculatePlanPositions();

        _searching = false;
        return true;
    }
    else
    {
        // target may be out range
        if (!InRange(xs, zs))
        {
#if LOG_POSITION_PROBL
            LOG_DEBUG(AI, "Problem: {}, start point {:.0f}, {:.0f} out of range", (const char*)_vehicle->GetDebugName(),
                      xs * LandGrid + 0.5 * LandGrid, zs * LandGrid + 0.5 * LandGrid);
#endif
            _searching = false;
            return false; // no solution
        }
        if (!InRange(xe, ze))
        {
#if LOG_POSITION_PROBL
            LOG_DEBUG(AI, "Problem: {}, end point {:.0f}, {:.0f} out of range", (const char*)_vehicle->GetDebugName(),
                      xe * LandGrid + 0.5 * LandGrid, ze * LandGrid + 0.5 * LandGrid);
#endif
            _searching = false;
            return false; // no solution
        }
        if (GetFieldCost(xs, zs) >= GET_UNACCESSIBLE)
        {
#if LOG_POSITION_PROBL
            LOG_DEBUG(AI, "Problem: {}, starting point {:.0f}, {:.0f} unaccessible",
                      (const char*)_vehicle->GetDebugName(), xs * LandGrid + 0.5 * LandGrid,
                      zs * LandGrid + 0.5 * LandGrid);
#endif
            if (!FindNearestSafe(xs, zs, 10000.0F))
            {
#if LOG_POSITION_PROBL
                LOG_DEBUG(AI, "Problem: {}, no accessible point found", (const char*)_vehicle->GetDebugName());
#endif
                _searching = false;
                return false; // no solution
            }
#if LOG_POSITION_PROBL
            LOG_DEBUG(AI, "Problem: {}, new starting point {:.0f}, {:.0f} accessible",
                      (const char*)_vehicle->GetDebugName(), xs * LandGrid + 0.5 * LandGrid,
                      zs * LandGrid + 0.5 * LandGrid);
#endif
        }
        if (GetFieldCost(xe, ze) >= GET_UNACCESSIBLE)
        {
#if LOG_POSITION_PROBL
            LOG_DEBUG(AI, "Problem: {}, end point {:.0f}, {:.0f} unaccessible", (const char*)_vehicle->GetDebugName(),
                      xe * LandGrid + 0.5 * LandGrid, ze * LandGrid + 0.5 * LandGrid);

#endif
            if (!FindNearestSafe(xe, ze, 10000.0F))
            {
#if LOG_POSITION_PROBL
                LOG_DEBUG(AI, "Problem: {}, no accessible point found", (const char*)_vehicle->GetDebugName());
#endif
                _searching = false;
                return false; // no solution
            }
#if LOG_POSITION_PROBL
            LOG_DEBUG(AI, "Problem: {}, end point {:.0f}, {:.0f} accessible", (const char*)_vehicle->GetDebugName(),
                      xe * LandGrid + 0.5 * LandGrid, ze * LandGrid + 0.5 * LandGrid);
#endif
            // new destination
            _destination[0] = xe * LandGrid + 0.5 * LandGrid;
            _destination[2] = ze * LandGrid + 0.5 * LandGrid;
            _destination[1] = 0; // not used
        }

        _heuristic = CalculateHeuristic(xs, zs, xe, ze);
        ASSField start(xs, zs);
        ASSField end(xe, ze);
        ASSCostFunction costFunction(_costFunction, _param, _vehicle);
        ASSHeuristicFunction heuristicFunction(_heuristic);
        _algorithm = new AStarStrategic(start, end, costFunction, heuristicFunction, GET_UNACCESSIBLE);

        // continue witch searching
        _searching = true;
        _iterTotal = 0;
#if CHECK_PERFORMANCE
        _timeTotal = 0;
#endif
        return true;
    }
}

void AIPathPlanner::StopSearching()
{
    _searching = false;
    _algorithm = nullptr;
}

bool AIPathPlanner::ProcessSearching()
{
    if (!_vehicle)
    {
        _searching = false;
#if LOG_STRAT
        LOG_DEBUG(AI, "Strategic path not found - vehicle destroyed");
#endif
        _algorithm = nullptr;
        return false;
    }

#if CHECK_PERFORMANCE
    _perfStart = ReadTsc();
#endif

    _iterTotal += _algorithm->Process(ITER_PER_CYCLE);
#if CHECK_PERFORMANCE
    _timeTotal += int(ReadTsc() - _perfStart);
#endif
    if (_algorithm->IsDone() || _iterTotal >= MAX_ITER)
    {
        _searching = false;

        const ASSNode* last = nullptr;
        bool partial = false;
        if (_algorithm->IsFound())
        {
            last = _algorithm->GetLastNode();
        }
        else
        {
            // destination not reached - search space or iteration budget exhausted
            // fall back to the most promising explored field and let the unit
            // replan from there (see alternate goal handling in OperMap)
            const ASSNode* best = _algorithm->GetBestNode();
            if (best && best->_parent)
            {
                const ASSNode* root = best;
                while (root->_parent)
                {
                    root = root->_parent;
                }
                // accept only when it brings us nearer to the destination than the start
                if (best->_f - best->_g < root->_f - root->_g)
                {
                    last = best;
                    partial = true;
                }
            }
        }

        if (last)
        {
            int depth = 0;
            for (const ASSNode* cur = last; cur != nullptr; cur = cur->_parent)
            {
                depth++;
            }
            _plan.Resize(depth);
            int i = depth - 1;
            for (const ASSNode* cur = last; cur != nullptr; cur = cur->_parent)
            {
                FieldPassing& info = _plan[i--];
                info._x = cur->_field.coord.x;
                info._z = cur->_field.coord.z;
                info._mode = FieldPassing::Move;
                info._cost = cur->_g;
            }

            if (partial)
            {
                // plan ends short of the wanted position
                // aim there - CalculatePlanPositions ends the plan in _destination
                _destination[0] = last->_field.coord.x * LandGrid + 0.5 * LandGrid;
                _destination[2] = last->_field.coord.z * LandGrid + 0.5 * LandGrid;
                _destination[1] = 0; // not used
            }

#if LOG_STRAT
            LOG_DEBUG(AI, "Strategic path found{}: {}, length {}, cost {:.0f} (in {} steps):",
                      partial ? " (partial)" : "", (const char*)_vehicle->GetDebugName(), depth, last->_g, _iterTotal);
#endif

            CalculatePlanPositions();

            _algorithm = nullptr;
            return true;
        }
        else
        {
#if LOG_STRAT
            LOG_DEBUG(AI, "Strategic path not found: {} (in {} steps)",
                      (const char*)_vehicle->GetDebugName(), _iterTotal);
#endif
            _algorithm = nullptr;
            return false;
        }
    }
    else
    {
        return true; // continue with searching
    }

}

float AIPathPlanner::CalculateHeuristic(int xs, int zs, int xe, int ze)
{
    float dx = xe - xs;
    float dz = ze - zs;
    int incz = zs < ze ? 1 : -1;
    int incx = xs < xe ? 1 : -1;
    float minCost = FLT_MAX;
    float maxCost = -FLT_MAX;
    float sumCost = 0;
    int count = 0;

    if (fabs(dx) < fabs(dz))
    {
        float invabsdz = 1.0 / fabs(dz);
        dx *= invabsdz;
        float x = xs + 0.5;
        for (int z = zs; z != ze; z += incz)
        {
            x += dx;
            if (toIntFloor(x) != xs)
            {
                AI_ERROR(toIntFloor(x) == xs + incx);
                float cost1 = GetFieldCost(xs, z);
                float cost2 = GetFieldCost(xs + incx, z);
                if (cost1 >= GET_UNACCESSIBLE || cost2 >= GET_UNACCESSIBLE)
                {
                    xs += incx;
                    continue;
                }
                saturateMax(maxCost, cost1);
                saturateMin(minCost, cost2);
                sumCost += cost1 + cost2;
                count += 2;

                xs += incx;
            }
            else
            {
                float cost = GetFieldCost(xs, z);
                if (cost >= GET_UNACCESSIBLE)
                {
                    continue;
                }
                saturateMax(maxCost, cost);
                saturateMin(minCost, cost);
                sumCost += cost;
                count++;
            }
        }
    }
    else if (fabs(dx) > 0)
    {
        float invabsdx = 1.0 / fabs(dx);
        dz *= invabsdx;
        float z = zs + 0.5;
        for (int x = xs; x != xe; x += incx)
        {
            z += dz;
            if (toIntFloor(z) != zs)
            {
                AI_ERROR(toIntFloor(z) == zs + incz);
                float cost1 = GetFieldCost(x, zs);
                float cost2 = GetFieldCost(x, zs + incz);
                if (cost1 >= GET_UNACCESSIBLE || cost2 >= GET_UNACCESSIBLE)
                {
                    zs += incz;
                    continue;
                }
                saturateMax(maxCost, cost1);
                saturateMin(minCost, cost2);
                sumCost += cost1 + cost2;
                count += 2;

                zs += incz;
            }
            else
            {
                float cost = GetFieldCost(x, zs);
                if (cost >= GET_UNACCESSIBLE)
                {
                    continue;
                }
                saturateMax(maxCost, cost);
                saturateMin(minCost, cost);
                sumCost += cost;
                count++;
            }
        }
    }
    else
    {
        // LOG_DEBUG(AI, "Singular search");
    }

    if (count <= 0)
    {
        maxCost = _vehicle->GetType()->GetMinCost();
        AI_ERROR(maxCost < GET_UNACCESSIBLE);
        maxCost *= LandGrid;
        minCost = sumCost = maxCost;
        count = 1;
    }
    return 0.9 * minCost + 0.1 * maxCost;
}

bool AIPathPlanner::FindNearestSafe(int& x, int& z, float threshold)
{
    if (IsSafe(x, z, threshold))
    {
        return true;
    }
    int rMax = 10; // do not too far (max 500 m)
    int i, r, xt, zt;
    for (r = 1; r < rMax; r++)
    {
        for (i = 0; i < r; i++)
        {
            xt = x - i;
            zt = z - r;
            if (IsSafe(xt, zt, threshold))
            {
                goto SafeFieldFound;
            }
            xt = x + i;
            zt = z - r;
            if (IsSafe(xt, zt, threshold))
            {
                goto SafeFieldFound;
            }
            xt = x + r;
            zt = z - i;
            if (IsSafe(xt, zt, threshold))
            {
                goto SafeFieldFound;
            }
            xt = x + r;
            zt = z + i;
            if (IsSafe(xt, zt, threshold))
            {
                goto SafeFieldFound;
            }
            xt = x + i;
            zt = z + r;
            if (IsSafe(xt, zt, threshold))
            {
                goto SafeFieldFound;
            }
            xt = x - i;
            zt = z + r;
            if (IsSafe(xt, zt, threshold))
            {
                goto SafeFieldFound;
            }
            xt = x - r;
            zt = z + i;
            if (IsSafe(xt, zt, threshold))
            {
                goto SafeFieldFound;
            }
            xt = x - r;
            zt = z - i;
            if (IsSafe(xt, zt, threshold))
            {
                goto SafeFieldFound;
            }
        }
        xt = x - r;
        zt = z - r;
        if (IsSafe(xt, zt, threshold))
        {
            goto SafeFieldFound;
        }
        xt = x + r;
        zt = z - r;
        if (IsSafe(xt, zt, threshold))
        {
            goto SafeFieldFound;
        }
        xt = x + r;
        zt = z + r;
        if (IsSafe(xt, zt, threshold))
        {
            goto SafeFieldFound;
        }
        xt = x - r;
        zt = z + r;
        if (IsSafe(xt, zt, threshold))
        {
            goto SafeFieldFound;
        }
    }
    return false;
SafeFieldFound:
    x = xt;
    z = zt;
    return true;
}

float AIPathPlanner::GetExposure(int x, int z)
{
    AI_ERROR(_vehicle);
    AIUnit* brain = _vehicle->CommanderUnit();
    if (!brain)
    {
        brain = _vehicle->PilotUnit();
    }
    AI_ERROR(brain);
    if (brain->IsHoldingFire())
    {
        return brain->GetGroup()->GetCenter()->GetExposurePessimistic(x, z);
    }
    else
    {
        return brain->GetGroup()->GetCenter()->GetExposureOptimistic(x, z);
    }
}

float AIPathPlanner::Distance(int index, Vector3Par pos) const
{
    Vector3 beg, end;
    GetPlanPosition(index, beg);
    GetPlanPosition(index + 1, end);
    Vector3 e = end - beg;
    Vector3 p = pos - beg;
    float tNext = (e * p) / e.SquareSizeXZ();
    saturate(tNext, 0, 1);
    Vector3 nearest = beg + tNext * e;
    float dist2Next = (nearest - pos).SquareSizeXZ();
    return dist2Next;
}

int AIPathPlanner::FindBestIndex(Vector3Par pos) const
{
    int n = _plan.Size();
    if (n == 0)
    {
        return -1;
    }

    int planIndex = -1;
    float dist2Min = FLT_MAX;
    for (int i = 0; i < n - 1; i++)
    {
        float dist2 = Distance(i, pos);
        if (dist2 < dist2Min)
        {
            dist2Min = dist2;
            planIndex = i;
        }
    }
    planIndex += 2;
    saturate(planIndex, 0, _plan.Size() - 1);
    return planIndex;
}

void AIPathPlanner::CalculatePlanPositions()
{
    int n = _plan.Size();
    _planPoints.Resize(n);
    if (n <= 0)
    {
        return;
    }

    VehicleWithAI* veh = _vehicle;
    AIUnit* unit = veh->CommanderUnit();
    if (!unit)
    {
        unit = veh->PilotUnit();
    }
    CombatMode cmode = unit ? unit->GetCombatMode() : CMAware;

    bool useRoads = cmode <= CMAware;

    for (int i = 0; i < n - 1; i++)
    {
        FieldPassing& pass = _plan[i];
        pass._mode = FieldPassing::Move;
        Vector3 pos;
        pos.Init();
        pos[0] = LandGrid * pass._x + 0.5 * LandGrid;
        pos[1] = 0;
        pos[2] = LandGrid * pass._z + 0.5 * LandGrid;
        if (useRoads)
        {
            GeographyInfo info = GLOB_LAND->GetGeography(pass._x, pass._z);
            if (info.u.road)
            { // try to find point on road
                float minDist2 = Square(40.0f);
                float found = false;
                Vector3 bestPos = pos;
                for (int zi = pass._z - 1; zi <= pass._z + 1; zi++)
                {
                    for (int xi = pass._x - 1; xi <= pass._x + 1; xi++)
                    {
                        if (!InRange(xi, zi))
                        {
                            continue;
                        }
                        RoadList& roadList = GRoadNet->GetRoadList(xi, zi);
                        for (int i = 0; i < roadList.Size(); i++)
                        {
                            RoadLink* item = roadList[i];
                            Vector3 pt = item->GetCenter();
                            float dist2 = (pos - pt).SquareSizeXZ();
                            if (dist2 < minDist2)
                            {
                                found = true;
                                bestPos = pt;
                                minDist2 = dist2;
                            }
                        }
                    }
                }
                if (found)
                {
                    pos = bestPos;
                    pass._mode = FieldPassing::MoveOnRoad;
                }
            }
        }
        pos[1] = GLOB_LAND->RoadSurfaceY(pos[0], pos[2]);
        _planPoints[i] = pos;
    }

    _planPoints[n - 1] = _destination;
    if (GRoadNet->IsOnRoad(_destination, 0))
    {
        _plan[n - 1]._mode = FieldPassing::MoveOnRoad;
    }
    else
    {
        _plan[n - 1]._mode = FieldPassing::Move;
    }
}

bool AIPathPlanner::GetPlanPosition(int index, Vector3& pos) const
{
    int n = _planPoints.Size();
    AI_ERROR(n == _plan.Size()); // _planPoints valid

    if (n <= 0)
    {
        pos = _destination;
        return true;
    }

    saturate(index, 0, n - 1);

    pos = _planPoints[index];
    return index == n - 1;
}

FieldPassing::Mode AIPathPlanner::GetPlanMode(int index) const
{
    if (index < 0 || index >= _plan.Size())
    {
        Fail("Out of plan");
        return FieldPassing::Move;
    }

    return (FieldPassing::Mode)_plan[index]._mode;
}

float AIPathPlanner::GetTotalCost() const
{
    int index = _plan.Size() - 1;
    if (index < 0)
    {
        return 0;
    }

    return _plan[index]._cost;
}

GeographyInfo AIPathPlanner::GetGeography(int index) const
{
    AI_ERROR(index >= 0 && index < _plan.Size());

    int x = _plan[index]._x;
    int z = _plan[index]._z;
    return GLOB_LAND->GetGeography(x, z);
}

bool AIPathPlanner::IsOnPath(int x, int z, int from, int to) const
{
    if (from < 0 || from >= _plan.Size())
    {
        return false;
    }
    if (to < 0 || to >= _plan.Size())
    {
        return false;
    }

    for (int i = from; i <= to; i++)
    {
        const FieldPassing& info = _plan[i];
        if (info._x == x && info._z == z)
        {
            return true;
        }
    }
    return false;
}

} // namespace Poseidon
