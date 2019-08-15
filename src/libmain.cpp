
#include "MatchVampPlugin.h"
#include "MatchTipicVampPlugin.h"

#include <vamp-sdk/PluginAdapter.h>

static Vamp::PluginAdapter<MatchVampPlugin> mvpAdapter;
static Vamp::PluginAdapter<MatchTipicVampPlugin> mtvpAdapter;

const VampPluginDescriptor *vampGetPluginDescriptor(unsigned int version,
                                                    unsigned int index)
{
    if (version < 1) return 0;

    switch (index) {
    case  0: return mvpAdapter.getDescriptor();
    case  1: return mtvpAdapter.getDescriptor();
    default: return 0;
    }
}
